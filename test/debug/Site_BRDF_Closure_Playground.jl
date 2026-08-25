# ============================== SETTINGS ==============================
# --- what to test ---
const WL          = 757.0     # wavelength (nm); SIF needs 640–850. Use 752 with CONSERVING = true.
const CI_FLAG     = true      # FEATURES.CI_IN_EXTINCTION_ONLY  (false = legacy ci prefactors on p_sensor / pso)
const STEM_RS     = false     # add the external stem-rescatter SIF source term (dob_stem/dof_stem); only matters if SAI > 0
const CONSERVING  = false     # zero leaf absorption in [SAC_LO,SAC_HI] + soil ρ = 1
const SAC_LO, SAC_HI = 750.0, 755.0
const HOTSPOT_Q   = nothing   # nothing = Emerald default (leaf width / canopy height ≈ 0.007); number to override (SCOPE: 0.05)
const N_MU, N_PHI = 24, 24    # hemisphere quadrature (24x24 agrees with 72x72 to <0.01%)
# --- canopy structure ---
const LAI         = 3.0
const CI          = 0.9
const SAI         = 0.0
const SZA         = 30.0
const LIDF_AB     = (0.0, 0.0)    # Verhoef (A,B): uniform (0,0) | planophile (1,0) | erectophile (-1,0) | plagiophile (0,-1) | extremophile (0,1) | spherical (-0.35,-0.15)
const N_LAYER_PER_LAI = 10        # layers = max(10, ceil(LAI*10))
# --- surfaces ---
const SOIL_RHO    = 0.15          # scalar shortwave soil albedo (prescribed)
const STEM_RHO    = (0.21, 0.49)  # stem reflectance VIS (<700 nm) / NIR (≥700 nm)
# --- leaf / physiology (forest-default baseline) ---
const CAB         = 35.0
const LMA         = 0.008
const MESO_N      = 1.4
const VCMAX       = 60.0
const JMAX_RATIO  = 1.73          # Jmax25 = JMAX_RATIO * Vcmax25
const VCMAX_EXPO  = 0.3
const G1          = 130.0         # Medlyn slope
const SPECTRA_SET = :new          # :new = NEW_PHI_2021, :old = OLD_PHI_2021
# ====================================================================

# cd("/home/klliu/SIFEsc/EmeraldTestV2"); using Pkg; Pkg.activate("."; io = devnull)
using Emerald
import Emerald.SPAC as ESPAC
import Emerald.Namespace as ENS
import Emerald.CanopyOptics as ECO
import Emerald.Land as ELND
const FT = Float64

# ── canopy builder (mirrors the forest-default spac_traits) ──────────────────────
function build()
    config = ENS.SPACConfig(FT; dataset = SPECTRA_SET == :new ? ENS.NEW_PHI_2021 : ENS.OLD_PHI_2021)
    #config.FEATURES.ENABLE_LEAF_SIF_SIGMOID = true
    #config.FEATURES.ENABLE_LEAF_SIF_RESCALE = true
    #config.FEATURES.ENABLE_CHL_SIF_SIGMOID = true
    #config.FEATURES.ENABLE_CHL_SIF_RESCALE = true
    config.METHODS.SOIL_ALBEDO = ENS.SoilAlbedoPrescribe()
    config.METHODS.STOMATAL_MODEL = ENS.MedlynSM{FT}()
    config.METHODS.STOMATAL_MODEL.G1 = G1
    config.FEATURES.EFFECTIVE_LEAF_SPECTRA = false
    #config.FEATURES.CI_IN_EXTINCTION_ONLY = CI_FLAG

    # conserving test: zero every pigment/water/dry-matter absorption cross-section in the band
    if CONSERVING
        Λ = config.CONSTANTS.SPECTRA.Λ; mask = (Λ .>= SAC_LO) .& (Λ .<= SAC_HI)
        for k in (:K_ANT, :K_BROWN, :K_CAB, :K_CAR_V, :K_CAR_Z, :K_CBC, :K_H₂O, :K_LMA, :K_PRO)
            getproperty(config.CONSTANTS.SPECTRA, k)[mask] .= 1e-6
        end
    end
    # stem reflectance
    Λ = config.CONSTANTS.SPECTRA.Λ
    config.CONSTANTS.SPECTRA.ρ_STEM[Λ .< 700] .= STEM_RHO[1]
    config.CONSTANTS.SPECTRA.ρ_STEM[Λ .>= 700] .= STEM_RHO[2]

    n_layer = max(N_LAYER_PER_LAI, Int(ceil(LAI * N_LAYER_PER_LAI)))
    spac = ENS.BulkSPAC(config; air_bounds = collect(0:(5.999 / (n_layer - 1)):13))
    ESPAC.prescribe_traits!(config, spac; lai = LAI, ci = CI, vcmax = VCMAX, jmax = JMAX_RATIO * VCMAX,
                            vertical_expo = VCMAX_EXPO, cab = CAB, sai = SAI)
    for leaf in spac.plant.leaves
        leaf.bio.trait.lma = LMA; leaf.bio.trait.meso_n = MESO_N
        leaf.capacitor.trait.v_max = 10.0; leaf.capacitor.state.v_storage = 10.0
    end
    spac.canopy.structure.trait.lidf.A = LIDF_AB[1]; spac.canopy.structure.trait.lidf.B = LIDF_AB[2]
    spac.canopy.sun_geometry.state.sza = SZA
    ESPAC.initialize_spac!(config, spac)
    spac.soil_bulk.auxil.ρ_sw .= CONSERVING ? 1.0 : SOIL_RHO
    ESPAC.spac!(config, spac, FT(3600))          # full solve (radiation, photosynthesis, ϕ_F)
    # re-run the canopy RT chain so geometry/hotspot settings are applied consistently
    geom!(config, spac); ECO.sensor_geometry!(config, spac)
    ECO.shortwave_radiation!(config, spac); ECO.reflection_spectrum!(config, spac); ECO.fluorescence_spectrum!(config, spac)
    return spac, config
end
geom!(config, spac) = isnothing(HOTSPOT_Q) ? ECO.sensor_geometry_aux!(config, spac) : ECO.sensor_geometry_aux!(config, spac.canopy, FT(HOTSPOT_Q))

# ── hemisphere quadrature: μ-midpoint × uniform azimuth; ∫∫ L cosθ sinθ dθ dφ = Σ L·μ·Δμ·Δφ ──
function hemisphere_nodes(n_mu, n_phi)
    μ = [(i - 0.5) / n_mu for i in 1:n_mu]; vza = acosd.(μ); φ = [(j - 0.5) * 360 / n_phi for j in 1:n_phi]
    return μ, vza, φ, (1 / n_mu) * (2π / n_phi)
end

# ── directional evaluations ──────────────────────────────────────────────────────
function set_view!(config, spac, vza, vaa)
    spac.canopy.sensor_geometry.state.vza = FT(vza); spac.canopy.sensor_geometry.state.vaa = FT(vaa)
    geom!(config, spac); ECO.sensor_geometry!(config, spac)
end
function refl_L(config, spac, vza, vaa, iwl)
    set_view!(config, spac, vza, vaa); ECO.reflection_spectrum!(config, spac)
    return spac.canopy.sensor_geometry.auxil.e_sensor[iwl]
end
function sif_terms(config, spac, vza, vaa, isif)
    set_view!(config, spac, vza, vaa); ECO.fluorescence_spectrum!(config, spac)
    se = spac.canopy.sensor_geometry.auxil; su = spac.canopy.sun_geometry.auxil; cs = spac.canopy.structure
    st = 0.0
    if STEM_RS && cs.trait.sai > 0
        IΛ = config.CONSTANTS.SPECTRA.IΛ_SIF; n = length(cs.trait.δlai)
        st = sum((se.dob_stem[IΛ,:][isif,:] .* su.e_sifꜜ[isif,1:n] .+ se.dof_stem[IΛ,:][isif,:] .* su.e_sifꜛ[isif,1:n]) .* se.p_sensor .* cs.trait.δsai) * se.ko_stem / pi
    end
    return (se.sif_obs_sunlit[isif], se.sif_obs_shaded[isif], se.sif_obs_scattered[isif] + st, se.sif_obs_soil[isif])
end
function integrate(f, config, spac)
    μ, vza, φ, w = hemisphere_nodes(N_MU, N_PHI); acc = nothing
    for j in eachindex(φ), i in eachindex(μ)
        v = f(config, spac, vza[i], φ[j]); acc = isnothing(acc) ? v .* (μ[i] * w) : acc .+ v .* (μ[i] * w)
    end
    return acc
end

# ── run ──────────────────────────────────────────────────────────────────────────
spac, config = build()
Λ = config.CONSTANTS.SPECTRA.Λ; ΛS = config.CONSTANTS.SPECTRA.Λ_SIF
iwl = argmin(abs.(Λ .- WL)); isif = argmin(abs.(ΛS .- WL))
su = spac.canopy.sun_geometry.auxil; rad = spac.meteo.rad_sw; saa = spac.canopy.sun_geometry.state.saa
leaf = spac.plant.leaves[end]

println("\n=== Canopy: LAI $(LAI)  CI $(CI)  SAI $(SAI)  SZA $(SZA)  LIDF $(LIDF_AB)  soil $(CONSERVING ? 1.0 : SOIL_RHO)  stem ρ $(STEM_RHO)  q=$(isnothing(HOTSPOT_Q) ? "default" : HOTSPOT_Q)  CI_FLAG=$(CI_FLAG)  STEM_RS=$(STEM_RS)  CONSERVING=$(CONSERVING)  layers=$(length(spac.plant.leaves))")
println("    wavelength: $(round(Λ[iwl],digits=1)) nm (SW grid) / $(round(ΛS[isif],digits=1)) nm (SIF grid);  leaf ρ=$(round(leaf.bio.auxil.ρ_leaf[iwl],digits=4)) τ=$(round(leaf.bio.auxil.τ_leaf[iwl],digits=4)) ρ+τ=$(round(leaf.bio.auxil.ρ_leaf[iwl]+leaf.bio.auxil.τ_leaf[iwl],digits=4))")

E_in = rad.e_dir[iwl] + rad.e_dif[iwl]; E_up = su.e_difꜛ[iwl,1]
L_int = integrate((c,s,v,a)->refl_L(c,s,v,a,iwl), config, spac)
L_nad = refl_L(config, spac, 0.0, saa, iwl)
println("\nREFLECTANCE @$(round(Λ[iwl],digits=0)) nm")
println("    albedo  E↑/E_in                 = ", round(E_up/E_in, digits=5))
println("    closure ∫L cosθ dΩ / E↑         = ", round(L_int/E_up, digits=4))
println("    nadir BRF  π·L_nadir/E_in       = ", round(pi*L_nad/E_in, digits=4), "    (π·L_nadir/E↑ = ", round(pi*L_nad/E_up, digits=4), ")")

E_sif = su.e_sifꜛ[isif,1]
T = integrate((c,s,v,a)->sif_terms(c,s,v,a,isif), config, spac)
Tn = sif_terms(config, spac, 0.0, saa, isif); Ln = sum(Tn)
println("\nSIF @$(round(ΛS[isif],digits=0)) nm")
println("    E_hemi (TOC upward SIF)         = ", round(E_sif, digits=4))
println("    closure ∫L cosθ dΩ / E_hemi     = ", round(sum(T)/E_sif, digits=4), "    [sunlit ", round(T[1]/E_sif,digits=3), " + shaded ", round(T[2]/E_sif,digits=3), " + rescattered ", round(T[3]/E_sif,digits=3), " + soil ", round(T[4]/E_sif,digits=3), "]")
println("    nadir SIF radiance              = ", round(Ln, digits=4), "    π·nadir/E_hemi = ", round(pi*Ln/E_sif, digits=4))
println("    escape TOC/leaf                 = ", round(E_sif / (sum(su.e_sifꜛ_layer[isif,:]) + sum(su.e_sifꜜ_layer[isif,:])), digits=4))
println("\nflux-side check (unchanged by CI_FLAG / HOTSPOT_Q): PPAR=", round(ELND.PPAR(spac),digits=1), "  GPP=", round(ELND.GPP(spac),digits=3), "  ΦF=", round(ELND.ΦF(spac),digits=5))

# This file contains function to compute the sun geometry related parameters (apart from sensor geometry)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function sun_geometry_aux! to update the trait-dependent auxiliary variables for sun geometry (to call in step_preparations!)
#     2024-Mar-01: compute the layer shortwave scattering fractions based on the new theory
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: use sza dependent clumping index
#     2024-Sep-09: account for diffuse CI impact on ks
# Bug fixes
#     2024-Mar-06: ci impact on fraction from solar direction
#
#######################################################################################################################################################################################################
"""

    sun_geometry_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update sun geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sun_geometry_aux! end;

sun_geometry_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = sun_geometry_aux!(config, spac.canopy);

sun_geometry_aux!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT} =
    sun_geometry_aux!(config, can.structure.trait, can.structure.t_aux, can.sun_geometry.state, can.sun_geometry.s_aux);

sun_geometry_aux!(config::SPACConfiguration{FT}, trait::CanopyStructureTrait{FT}, t_aux::CanopyStructureTDAuxil{FT}, sunst::SunGeometryState{FT}, sunsa::SunGeometrySDAuxil{FT}) where {FT} = (
    # if sza > 89 or both LAI and SAI are zero, do nothing
    if sunst.sza > 89 || (trait.lai <= 0 && trait.sai <= 0)
        return nothing
    end;

    (; Θ_AZI, Θ_INCL) = config;

    # compute the clumping index from solar zenith angle
    sunsa.ci_sun = trait.ci.ci_0 * (1 - trait.ci.ci_1 * cosd(sunst.sza));

    # extinction coefficients for the solar radiation
    for i in eachindex(Θ_INCL)
        Cs = cosd(Θ_INCL[i]) * cosd(sunst.sza);
        Ss = sind(Θ_INCL[i]) * sind(sunst.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        sunsa.Cs_incl[i] = Cs;
        sunsa.Ss_incl[i] = Ss;
        sunsa.βs_incl[i] = βs;
        sunsa.ks_incl[i] = 2 / FT(π) / cosd(sunst.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;
    sunsa.ks_leaf = t_aux.p_incl_leaf' * sunsa.ks_incl * sunsa.ci_sun;
    sunsa.ks_stem = t_aux.p_incl_stem' * sunsa.ks_incl * sunsa.ci_sun;

    # compute the scattering weights for diffuse/direct -> diffuse for backward and forward scattering
    sunsa.sdb_leaf = 0;
    sunsa.sdf_leaf = 0;
    sunsa.sdb_stem = 0;
    sunsa.sdf_stem = 0;
    for i in eachindex(Θ_INCL)
        f_ada = f_adaxial(sunst.sza, Θ_INCL[i]);
        f_aba = 1 - f_ada;
        f_inc = Θ_INCL[i] / 180;
        sunsa.sdb_leaf += (f_ada * (1 - f_inc) + f_aba * f_inc) * t_aux.p_incl_leaf[i];
        sunsa.sdf_leaf += (f_ada * f_inc + f_aba * (1 - f_inc)) * t_aux.p_incl_leaf[i];
        sunsa.sdb_stem += (f_ada * (1 - f_inc) + f_aba * f_inc) * t_aux.p_incl_stem[i];
        sunsa.sdf_stem += (f_ada * f_inc + f_aba * (1 - f_inc)) * t_aux.p_incl_stem[i];
    end;

    # compute the sunlit leaf fraction
    # ps(x) = sunsa.ci_sun * exp.(sunsa.ks .* trait.lai .* t_aux.x_bnds);
    kscipai = sunsa.ks_leaf * trait.lai + sunsa.ks_stem * trait.sai;
    for i in eachindex(trait.δlai)
        ksciipai = sunsa.ks_leaf * trait.δlai[i] + sunsa.ks_stem * trait.δsai[i];
        sunsa.p_sunlit[i] = sunsa.ci_sun / ksciipai * (exp(kscipai * t_aux.x_bnds[i]) - exp(kscipai * t_aux.x_bnds[i+1]));
    end;

    # compute the fs and fs_abs matrices
    for i in eachindex(Θ_AZI)
        view(sunsa.fs,:,i) .= sunsa.Cs_incl .+ sunsa.Ss_incl .* cosd(Θ_AZI[i]);
    end;
    sunsa.fs ./= cosd(sunst.sza);
    sunsa.fs_abs .= abs.(sunsa.fs);
    mul!(sunsa.fs_abs_mean, sunsa.fs_abs', t_aux.p_incl_leaf);
    for i in eachindex(Θ_INCL)
        view(sunsa.fs_cos²_incl,i,:) .= view(sunsa.fs,i,:) .* (cosd(Θ_INCL[i]) ^ 2);
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sun_geometry! (run per solar zenith angle)
#     2023-Oct-11: compute canopy layer scattering, reflectance, and transmittance
#     2023-Oct-14: do nothing if sza > 89 or LAI <= 0
#     2023-Oct-18: account for SAI in the sun geometry calculation
#     2024-Mar-01: compute the layer shortwave scattering coefficients based on the new theory
#     2024-Jun-06: fix a typo in the calculation of ρ_sd
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: use sza dependent clumping index
#     2024-Oct-16: add option to compute effective leaf spectra based on CI
#     2024-Nov-08: when using EFFECTIVE_LEAF_SPECTRA make sure LAI > 0
#
#######################################################################################################################################################################################################
"""

    sun_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update sun geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sun_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;
    sun_geo = spac.canopy.sun_geometry;

    if sun_geo.state.sza > 89 || (can_str.trait.lai <= 0 && can_str.trait.sai <= 0)
        return nothing
    end;

    # if sza <= 89 and LAI+SAI > 0, run the sun geometry function
    (; EFFECTIVE_LEAF_SPECTRA, SPECTRA) = config;
    leaves = spac.plant.leaves;
    sbulk = spac.soil_bulk;
    n_layer = length(leaves);

    # compute the effective leaf reflectance and transmittance spectra using the PROSPECT model scheme
    mask_effective = EFFECTIVE_LEAF_SPECTRA && can_str.trait.lai > 0;
    if mask_effective
        ρ_2 = spac.cache.cache_wl_1;
        τ_2 = spac.cache.cache_wl_2;
        n_eff = 1 / sun_geo.s_aux.ci_sun;
        for irt in 1:n_layer
            ilf = n_layer + 1 - irt;
            leaf = leaves[ilf];
            ρ_1 = leaf.bio.auxil.ρ_leaf;
            τ_1 = leaf.bio.auxil.τ_leaf;
            ρ_2 .= layer_2_ρ.(ρ_1, τ_1, n_eff - 1);
            τ_2 .= layer_2_τ.(ρ_1, τ_1, n_eff - 1);
            sun_geo.auxil.ρ_leaf_eff[:,irt] .= leaf_ρ.(ρ_1, τ_1, ρ_1, τ_1, ρ_2);
            sun_geo.auxil.τ_leaf_eff[:,irt] .= leaf_τ.(τ_1, ρ_1, ρ_2, τ_2);
        end;
    end;

    # compute the scattering coefficients for the solar radiation per leaf area
    # use effective leaf spectra if EFFECTIVE_LEAF_SPECTRA is true
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        ρ_leaf = mask_effective ? view(sun_geo.auxil.ρ_leaf_eff,:,irt) : leaf.bio.auxil.ρ_leaf;
        τ_leaf = mask_effective ? view(sun_geo.auxil.τ_leaf_eff,:,irt) : leaf.bio.auxil.τ_leaf;
        sun_geo.auxil.sdb_leaf[:,irt] .= sun_geo.s_aux.sdb_leaf .* ρ_leaf .+ sun_geo.s_aux.sdf_leaf .* τ_leaf;
        sun_geo.auxil.sdf_leaf[:,irt] .= sun_geo.s_aux.sdf_leaf .* ρ_leaf .+ sun_geo.s_aux.sdb_leaf .* τ_leaf;
        sun_geo.auxil.sdb_stem[:,irt] .= sun_geo.s_aux.sdb_stem .* SPECTRA.ρ_STEM;
        sun_geo.auxil.sdf_stem[:,irt] .= sun_geo.s_aux.sdf_stem .* SPECTRA.ρ_STEM;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    # Similarly, we computed the transmittance and reflectance for the whole layer using expotential functions
    #     sun_geo.auxil.τ_ss_layer .= exp.(-1 .* sun_geo.s_aux.ks .* can_str.trait.δlai .* sun_geo.s_aux.ci_sun);
    #     sun_geo.auxil.τ_sd_layer .= 1 .- exp.(-1 .* sun_geo.auxil.sdf_leaf .* can_str.trait.δlai' .* sun_geo.s_aux.ci_sun);
    #     sun_geo.auxil.ρ_sd_layer .= 1 .- exp.(-1 .* sun_geo.auxil.sdb_leaf .* can_str.trait.δlai' .* sun_geo.s_aux.ci_sun);
    # Later, we included SAI as well.
    # However, as of 2024-Feb-29, we found an issue with the equations above when LAI and SAI are big enough in a single layer that sum of reflectance and transmittance is greater than 1.
    # Therefore, we revised the equations using calculus and the equations for the whole layer are as follows:
    #     sun_geo.auxil.ρ_sd_layer = ∫_0^iCIPAI (sdb_leaf * δLAI + sbd_stem * δSAI) / δPAI * ks * exp(-ks * x) * dx = (sdb_leaf * δLAI + sbd_stem * δSAI) / δPAI * (1 - exp(-ks * iCIPAI))
    #     sun_geo.auxil.τ_sd_layer = ∫_0^iCIPAI (sdf_leaf * δLAI + sdf_stem * δSAI) / δPAI * ks * exp(-ks * x) * dx = (sdf_leaf * δLAI + sdf_stem * δSAI) / δPAI * (1 - exp(-ks * iCIPAI))
    # Similarly, we used the same logic for sensor geometry and canopy structure.
    # As of 2014-Sep-09: we accounted for CI impact through ks_leaf and ks_stem, so we do not need to multiply ks by ci_sun in the code below.
    kt_sd_x = spac.cache.cache_wl_1;
    kr_sd_x = spac.cache.cache_wl_2;
    for i in 1:n_layer
        δlai = can_str.trait.δlai[i];
        δsai = can_str.trait.δsai[i];
        δpai = δlai + δsai;
        kt_ss_x = sun_geo.s_aux.ks_leaf * δlai + sun_geo.s_aux.ks_stem * δsai;
        kt_sd_x .= (view(sun_geo.auxil.sdf_leaf,:,i) .* δlai .+ view(sun_geo.auxil.sdf_stem,:,i) .* δsai) ./ δpai;
        kr_sd_x .= (view(sun_geo.auxil.sdb_leaf,:,i) .* δlai .+ view(sun_geo.auxil.sdb_stem,:,i) .* δsai) ./ δpai;
        sun_geo.auxil.τ_ss_layer[i] = exp(-kt_ss_x);
        sun_geo.auxil.τ_sd_layer[:,i] .= (1 - sun_geo.auxil.τ_ss_layer[i]) .* kt_sd_x;
        sun_geo.auxil.ρ_sd_layer[:,i] .= (1 - sun_geo.auxil.τ_ss_layer[i]) .* kr_sd_x;
    end;

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    sun_geo.auxil.ρ_sd[:,end] .= sbulk.auxil.ρ_sw;
    for i in n_layer:-1:1
        ρ_dd_layer = view(can_str.auxil.ρ_dd_layer,:,i  );
        ρ_dd_j     = view(can_str.auxil.ρ_dd      ,:,i+1);
        ρ_sd_layer = view(sun_geo.auxil.ρ_sd_layer,:,i  );
        ρ_sd_i     = view(sun_geo.auxil.ρ_sd      ,:,i  );
        ρ_sd_j     = view(sun_geo.auxil.ρ_sd      ,:,i+1);
        τ_dd_i     = view(can_str.auxil.τ_dd      ,:,i  );
        τ_sd_layer = view(sun_geo.auxil.τ_sd_layer,:,i  );
        τ_sd_i     = view(sun_geo.auxil.τ_sd      ,:,i  );
        τ_ss_layer = view(sun_geo.auxil.τ_ss_layer,  i  );

        τ_sd_i .= (τ_sd_layer .+ τ_ss_layer .* ρ_sd_j .* ρ_dd_layer) ./ (1 .- ρ_dd_layer .* ρ_dd_j);    # sdit + ssit-sdjr-ddit; rescale
        ρ_sd_i .= ρ_sd_layer .+ τ_ss_layer .* ρ_sd_j .* τ_dd_i .+ τ_sd_layer .* ρ_dd_j .* τ_dd_i;       # sdir + ssit-sdjr-ddit + sdit-ddjr-ddit
    end;

    return nothing
end;

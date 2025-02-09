module LeafOptics

using SpecialFunctions: expint

using ..EmeraldPhysics.Constant: M_H₂O, ρ_H₂O
using ..EmeraldUtility.StructEqual: sync_struct!

using ..Namespace: SIFMatrixDualspectMethod, SIFMatrixFluspectMethod, SIFMatrixPlatespectMethod
using ..Namespace: LeafBio, LeafBioState, LeafBioTrait
using ..Namespace: BulkSPAC, SPACConfiguration


include("prospect/interface.jl");
include("prospect/sublayer.jl");
include("prospect/layer.jl");
include("prospect/leaf.jl");

include("platespect/effective.jl");
include("platespect/fluorescence.jl");

# Both Fluspect and Dualspect use the same doubling adding method
include("kubelka-munk/fluorescence.jl");

include("dualspect/fluorescence.jl");
include("fluspect/fluorescence.jl");


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to run all the step within one function all
#     2023-Sep-16: compute SIF conversion matrices within this function
#     2024-Aug-13: add option to compute SIF matrices using N = nothing (default) for integral method (super fast and accurate) or N = Int for numerical method
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, lwc::FT, θ::FT = FT(40)) where {FT}

Update the interface, sublayer, layer, and leaf level reflectance and transmittance within `bio`, given
- `config` SPAC configuration
- `bio` LeafBio struct
- `lwc` Leaf water content
- `θ` Incoming radiation angle

"""
function leaf_spectra!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, lwc::FT, θ::FT = FT(40)) where {FT}
    leaf_interface_ρ_τ!(config, bio, θ);
    leaf_sublayer_f_τ!(config, bio, lwc);
    leaf_layer_ρ_τ!(bio);
    leaf_ρ_τ!(bio);

    if config.ENABLE_SIF
        leaf_sif_matrices!(config, bio);
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2022-Jun-29: add method for BulkSPAC
#     2024-Feb-22: do nothing if lai is zero
#     2024-Feb-28: support VERTICAL_BIO feature to update leaf spectra from top leaf
#
#######################################################################################################################################################################################################
"""

    plant_leaf_spectra!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update leaf reflectance and transmittance for SPAC, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function plant_leaf_spectra!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # if there is no leaf, do nothing
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    # update leaf reflectance and transmittance only if LAI > 0
    # use top leaf to update the rest
    (; VERTICAL_BIO) = config;
    if !VERTICAL_BIO
        topleaf = spac.plant.leaves[end];
        leaf_spectra!(config, topleaf.bio, topleaf.capacitor.state.v_storage);
        for i in 1:length(spac.plant.leaves)-1
            sync_struct!(topleaf.bio.auxil, spac.plant.leaves[i].bio.auxil);
        end;

        return nothing
    end;

    # update all leaves
    for leaf in spac.plant.leaves
        leaf_spectra!(config, leaf.bio, leaf.capacitor.state.v_storage);
    end;

    return nothing
end;


#=

#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: add function to compute leaf level PAR and APAR
#     2022-Jun-27: refactor the function to return PAR, APAR, and PPAR
#
#######################################################################################################################################################################################################
"""

    leaf_PAR(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::ShortwaveRadiation{FT}; apar_car::Bool = true) where {FT}

Return leaf level PAR, APAR, and PPAR, given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `spectra` `ReferenceSpectra` type struct that contains absorption characteristic curves
- `rad` `ShortwaveRadiation` type struct that contains incoming radiation information
- `apar_car` If true (default), account carotenoid absorption as PPAR; otherwise, PPAR is only by chlorophyll

---
# Examples
```julia
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
spectra = EmeraldNamespace.ReferenceSpectra{Float64}();
rad = EmeraldNamespace.ShortwaveRadiation{Float64}();
par,apar,ppar = leaf_PAR(bio, spectra, rad);
par,apar,ppar = leaf_PAR(bio, spectra, rad; apar_car=false);
```

"""
function leaf_PAR(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::ShortwaveRadiation{FT}; apar_car::Bool = true) where {FT}
    (; IΛ_PAR, ΔΛ_PAR, Λ_PAR) = spectra;

    # PPAR absorption feature (after APAR is computed)
    _α_ppar = (apar_car ? view(bio.α_cabcar, IΛ_PAR) : view(bio.α_cab, IΛ_PAR));

    # PAR, APAR, and PPAR energy from direct and diffuse light
    _e_par_dir   = view(rad.e_dir, IΛ_PAR);
    _e_par_diff  = view(rad.e_dif, IΛ_PAR);
    _e_apar_dir  = view(bio.α_sw, IΛ_PAR) .* _e_par_dir;
    _e_apar_diff = view(bio.α_sw, IΛ_PAR) .* _e_par_diff;
    _e_ppar_dir  = _α_ppar .* _e_apar_dir;
    _e_ppar_diff = _α_ppar .* _e_apar_diff;

    # PAR, APAR, and PPAR photons from direct and diffuse light
    _par_dir   = photon.(Λ_PAR, _e_par_dir  );
    _par_diff  = photon.(Λ_PAR, _e_par_diff );
    _apar_dir  = photon.(Λ_PAR, _e_apar_dir );
    _apar_diff = photon.(Λ_PAR, _e_apar_diff);
    _ppar_dir  = photon.(Λ_PAR, _e_ppar_dir );
    _ppar_diff = photon.(Λ_PAR, _e_ppar_diff);

    # total PAR and APAR in μmol photons m⁻² s⁻¹
    _Σpar_dir   = _par_dir'   * ΔΛ_PAR * 1000;
    _Σpar_diff  = _par_diff'  * ΔΛ_PAR * 1000;
    _Σapar_dir  = _apar_dir'  * ΔΛ_PAR * 1000;
    _Σapar_diff = _apar_diff' * ΔΛ_PAR * 1000;
    _Σppar_dir  = _ppar_dir'  * ΔΛ_PAR * 1000;
    _Σppar_diff = _ppar_diff' * ΔΛ_PAR * 1000;

    return _Σpar_dir + _Σpar_diff, _Σapar_dir + _Σapar_diff, _Σppar_dir + _Σppar_diff
end;


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Jul-08: add leaf level SIF simulation
#     2021-Jul-08: use mat_b and mat_f for SIF at backward and forward directions
#     2021-Aug-05: add option to sumulate SIF in photon to photon mode
#     2021-Oct-22: refactor the function to leaf_SIF to return the SIFs directly
#
#######################################################################################################################################################################################################
"""

    leaf_SIF(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::ShortwaveRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT}

Return the leaf level SIF at backward and forward directions, given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `spectra` `ReferenceSpectra` type struct that contains absorption characteristic curves
- `rad` `ShortwaveRadiation` type struct that contains incoming radiation information
- `ϕ` Fluorescence quantum yield
- `ϕ_photon` If true (default), convert photon to photon when computing SIF; otherwise, convert energy to energy

---
# Examples
```julia
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
spectra = EmeraldNamespace.ReferenceSpectra{Float64}();
rad = EmeraldNamespace.ShortwaveRadiation{Float64}();
sif_b,sif_f = leaf_SIF(bio, spectra, 0.01);
sif_b,sif_f = leaf_SIF(bio, spectra, 0.01; ϕ_photon=false);
```

"""
function leaf_SIF(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::ShortwaveRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT}
    (; IΛ_SIFE, ΔΛ_SIFE, Λ_SIF, Λ_SIFE) = spectra;

    # calculate the excitation energy and photons
    _e_excitation = (view(rad.e_dir, IΛ_SIFE) .+ view(rad.e_dif, IΛ_SIFE)) .* ΔΛ_SIFE;

    # convert energy to energy using the matrices
    if !ϕ_photon
        _sif_b = bio.mat_b * _e_excitation * ϕ / FT(π);
        _sif_f = bio.mat_f * _e_excitation * ϕ / FT(π);

        return _sif_b, _sif_f
    end;

    # convert energy to photon
    _phot_excitation = photon.(Λ_SIFE, _e_excitation);

    # convert photon to photon using the matrices
    _phot_b = bio.mat_b * _phot_excitation * ϕ / FT(π);
    _phot_f = bio.mat_f * _phot_excitation * ϕ / FT(π);

    # convert photon to back to energy
    _sif_b = energy.(Λ_SIF, _phot_b);
    _sif_f = energy.(Λ_SIF, _phot_f);

    return _sif_b, _sif_f
end;

=#


end; # module

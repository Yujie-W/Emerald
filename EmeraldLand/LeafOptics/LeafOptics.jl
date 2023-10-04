module LeafOptics

using SpecialFunctions: expint

using ..EmeraldPhysics.Constant: M_H₂O, ρ_H₂O
using ..EmeraldPhysics.Optics: energy, energy!, photon, photon!

using ..Namespace: HyperspectralRadiation, ReferenceSpectra
using ..Namespace: MultiLayerSPAC, SPACConfiguration


using ..Namespace: HyperLeafBio, HyperLeafBioState
include("doubling.jl");
include("effective.jl");
include("fluorescence.jl");
include("interface.jl");
include("layer.jl");
include("leaf.jl");
include("sublayer.jl");


#=
#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Oct-22: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT}

Update leaf reflectance and transmittance (e.g., prescribe broadband PAR and NIR values), given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `spectra` `ReferenceSpectra` type struct that contains absorption characteristic curves
- `ρ_par` Reflectance at PAR region
- `ρ_nir` Reflectance at NIR region
- `τ_par` Transmittance at PAR region
- `τ_nir` Transmittance at NIR region

# Examples
```julia
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
spectra = EmeraldNamespace.ReferenceSpectra{Float64}();
leaf_spectra!(bio, spectra, 0.1, 0.45, 0.05, 0.25);
```

"""
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT} = (
    (; IΛ_NIR, IΛ_PAR) = spectra;

    bio.ρ_sw[IΛ_PAR] .= ρ_par;
    bio.ρ_sw[IΛ_NIR] .= ρ_nir;
    bio.τ_sw[IΛ_PAR] .= τ_par;
    bio.τ_sw[IΛ_NIR] .= τ_nir;

    bio.α_sw = 1 .- bio.τ_sw .- bio.ρ_sw;

    return nothing
);
=#


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2022-Jun-29: add method for MultiLayerSPAC
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Update leaf reflectance and transmittance for SPAC, given
- `config` Configurations of spac model
- `spac` `MultiLayerSPAC` type SPAC

"""
leaf_spectra!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; LEAVES) = spac;

    for _leaf in LEAVES
        leaf_spectra!(config, _leaf.bio, _leaf.capacitor.state.v_storage);
    end;

    return nothing
);


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

    leaf_PAR(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::HyperspectralRadiation{FT}; apar_car::Bool = true) where {FT}

Return leaf level PAR, APAR, and PPAR, given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `spectra` `ReferenceSpectra` type struct that contains absorption characteristic curves
- `rad` `HyperspectralRadiation` type struct that contains incoming radiation information
- `apar_car` If true (default), account carotenoid absorption as PPAR; otherwise, PPAR is only by chlorophyll

---
# Examples
```julia
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
spectra = EmeraldNamespace.ReferenceSpectra{Float64}();
rad = EmeraldNamespace.HyperspectralRadiation{Float64}();
par,apar,ppar = leaf_PAR(bio, spectra, rad);
par,apar,ppar = leaf_PAR(bio, spectra, rad; apar_car=false);
```

"""
function leaf_PAR(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::HyperspectralRadiation{FT}; apar_car::Bool = true) where {FT}
    (; IΛ_PAR, ΔΛ_PAR, Λ_PAR) = spectra;

    # PPAR absorption feature (after APAR is computed)
    _α_ppar = (apar_car ? view(bio.α_cabcar, IΛ_PAR) : view(bio.α_cab, IΛ_PAR));

    # PAR, APAR, and PPAR energy from direct and diffuse light
    _e_par_dir   = view(rad.e_direct , IΛ_PAR);
    _e_par_diff  = view(rad.e_diffuse, IΛ_PAR);
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
end


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

    leaf_SIF(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::HyperspectralRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT}

Return the leaf level SIF at backward and forward directions, given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `spectra` `ReferenceSpectra` type struct that contains absorption characteristic curves
- `rad` `HyperspectralRadiation` type struct that contains incoming radiation information
- `ϕ` Fluorescence quantum yield
- `ϕ_photon` If true (default), convert photon to photon when computing SIF; otherwise, convert energy to energy

---
# Examples
```julia
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
spectra = EmeraldNamespace.ReferenceSpectra{Float64}();
rad = EmeraldNamespace.HyperspectralRadiation{Float64}();
sif_b,sif_f = leaf_SIF(bio, spectra, 0.01);
sif_b,sif_f = leaf_SIF(bio, spectra, 0.01; ϕ_photon=false);
```

"""
function leaf_SIF(bio::HyperspectralLeafBiophysics{FT}, spectra::ReferenceSpectra{FT}, rad::HyperspectralRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT}
    (; IΛ_SIFE, ΔΛ_SIFE, Λ_SIF, Λ_SIFE) = spectra;

    # calculate the excitation energy and photons
    _e_excitation = (view(rad.e_direct, IΛ_SIFE) .+ view(rad.e_diffuse, IΛ_SIFE)) .* ΔΛ_SIFE;

    # convert energy to energy using the matrices
    if !ϕ_photon
        _sif_b = bio.mat_b * _e_excitation * ϕ / FT(pi);
        _sif_f = bio.mat_f * _e_excitation * ϕ / FT(pi);

        return _sif_b, _sif_f
    end;

    # convert energy to photon
    _phot_excitation = photon.(Λ_SIFE, _e_excitation);

    # convert photon to photon using the matrices
    _phot_b = bio.mat_b * _phot_excitation * ϕ / FT(pi);
    _phot_f = bio.mat_f * _phot_excitation * ϕ / FT(pi);

    # convert photon to back to energy
    _sif_b = energy.(Λ_SIF, _phot_b);
    _sif_f = energy.(Λ_SIF, _phot_f);

    return _sif_b, _sif_f
end;

=#


end # module

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-Mar-13: add function to compute leaf level broadband reflectance and transmittance
#     2025-Mar-13: add function to compute canopy level albedo
#
#######################################################################################################################################################################################################
"""

    LEAF_RT_BB(config::SPACConfiguration{FT}, lbio::LeafBio{FT}) where {FT}

Return leaf PAR and NIR reflectance and transmittance from leaf level spectra, given
- `config`: a `SPACConfiguration` object containing the reference radiation spectra
- `lbio`: a `LeafBio` object containing the leaf level reflectance and transmittance spectra

"""
function LEAF_RT_BB end;

LEAF_RT_BB(config::SPACConfiguration{FT}, lbio::LeafBio{FT}) where {FT} = (
    (; SPECTRA) = config;

    # compute the par and nir radiation (reference radiation spectra)
    sunrad = SPECTRA.SOLAR_RAD[:,1] .+ SPECTRA.SOLAR_RAD[:,2];
    parrad = sunrad[SPECTRA.IΛ_PAR]' * SPECTRA.ΔΛ[SPECTRA.IΛ_PAR];
    nirrad = sunrad[SPECTRA.IΛ_NIR]' * SPECTRA.ΔΛ[SPECTRA.IΛ_NIR];

    # compute the leaf level reflectance
    par_refl = (lbio.auxil.ρ_leaf[SPECTRA.IΛ_PAR] .* sunrad[SPECTRA.IΛ_PAR])' * SPECTRA.ΔΛ[SPECTRA.IΛ_PAR] / parrad;
    par_tran = (lbio.auxil.τ_leaf[SPECTRA.IΛ_PAR] .* sunrad[SPECTRA.IΛ_PAR])' * SPECTRA.ΔΛ[SPECTRA.IΛ_PAR] / parrad;
    nir_refl = (lbio.auxil.ρ_leaf[SPECTRA.IΛ_NIR] .* sunrad[SPECTRA.IΛ_NIR])' * SPECTRA.ΔΛ[SPECTRA.IΛ_NIR] / nirrad;
    nir_tran = (lbio.auxil.τ_leaf[SPECTRA.IΛ_NIR] .* sunrad[SPECTRA.IΛ_NIR])' * SPECTRA.ΔΛ[SPECTRA.IΛ_NIR] / nirrad;

    return par_refl, par_tran, nir_refl, nir_tran
);

"""

    ALBEDO_BB(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the canopy level broadband albedo from the SPAC model, given
- `config`: a `SPACConfiguration` object containing the reference radiation spectra
- `spac`: a `BulkSPAC` object containing the canopy level albedo and the actual meteo radiation

"""
function ALBEDO_BB end;

ALBEDO_BB(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;

    # compute the par and nir radiation (using the actral meteo radiation)
    sunrad = spac.meteo.rad_sw.e_dif .+ spac.meteo.rad_sw.e_dir;
    parrad = sunrad[SPECTRA.IΛ_PAR]' * SPECTRA.ΔΛ[SPECTRA.IΛ_PAR];
    nirrad = sunrad[SPECTRA.IΛ_NIR]' * SPECTRA.ΔΛ[SPECTRA.IΛ_NIR];

    # compute the canopy level reflectance and transmittance
    par_albedo = (spac.canopy.sun_geometry.auxil.albedo[SPECTRA.IΛ_PAR] .* sunrad[SPECTRA.IΛ_PAR])' * SPECTRA.ΔΛ[SPECTRA.IΛ_PAR] / parrad;
    nir_albedo = (spac.canopy.sun_geometry.auxil.albedo[SPECTRA.IΛ_NIR] .* sunrad[SPECTRA.IΛ_NIR])' * SPECTRA.ΔΛ[SPECTRA.IΛ_NIR] / nirrad;

    return par_albedo, nir_albedo
);

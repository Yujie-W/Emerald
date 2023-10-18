#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-13: add function to add up total SIF photons (unit mol m⁻² s⁻¹)
#
#######################################################################################################################################################################################################
"""

ΣSIF(spac::BulkSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF end;

ΣSIF(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in photons unit
    Σsif::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        Σsif += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.ppar_sunlit .* leaves[i].flux.auxil.ϕ_f_sunlit) +
                (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.ppar_shaded * leaves[i].flux.auxil.ϕ_f_shaded) * canopy.structure.state.δlai[j];
    end;

    return Σsif
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-11: add function to add up total SIF at chloroplast level (without any reabsorption)
#     2023-Sep-11: convert the unit to W m⁻² per ground area from nW m⁻²
#
#######################################################################################################################################################################################################
"""

    ΣSIF_CHL(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF_CHL end;

ΣSIF_CHL(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in energy unit before reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(leaves)
        Σsif += view(canopy.sun_geometry.auxil.e_sif_chl,:,i)' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-11: add function to add up total SIF in energy unit at leaf level after reabsorption
#     2023-Sep-11: convert the unit to W m⁻² per ground area from nW m⁻²
#
#######################################################################################################################################################################################################
"""

    ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at leaf level after reabsorption in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF_LEAF end;

ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in energy unit after reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(leaves)
        Σsif += (canopy.sun_geometry.auxil.e_sifꜜ_layer[:,i] .+ canopy.sun_geometry.auxil.e_sifꜛ_layer[:,i])' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);

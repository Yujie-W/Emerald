# This file contains functions to compute the sensor geometry of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sensor_geometry! (run per viewing zenith angle)
#     2023-Oct-11: compute canopy layer scattering
#     2023-Oct-13: improve p_sun_sensor calculation accuracy
#     2023-Oct-14: if none of REF or SIF is enabled, skip the sensor geometry calculation
#     2023-Oct-18: account for SAI in the sensor geometry calculation
#     2024-Feb-22: add solar zenith angle control
#     2024-Feb-25: move the trait- and state-dependent calculations to the s_aux! function
#
#######################################################################################################################################################################################################
"""

    sensor_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sensor_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;
    sun_geo = spac.canopy.sun_geometry;

    if (!config.ENABLE_REF && !config.ENABLE_SIF) || sun_geo.state.sza > 89 || (can_str.trait.lai <= 0 && can_str.trait.sai <= 0)
        return nothing
    end;

    # run the sensor geometry simulations only if any of canopy reflectance feature or fluorescence feature is enabled and if LAI+SAI > 0
    (; SPECTRA) = config;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    n_layer = length(leaves);

    # compute the scattering coefficients per leaf area
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sen_geo.auxil.dob_leaf[:,irt] .= sen_geo.s_aux.dob * leaf.bio.auxil.ρ_leaf .+ sen_geo.s_aux.dof * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.dof_leaf[:,irt] .= sen_geo.s_aux.dof * leaf.bio.auxil.ρ_leaf .+ sen_geo.s_aux.dob * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.so_leaf[:,irt]  .= sen_geo.s_aux.sob * leaf.bio.auxil.ρ_leaf .+ sen_geo.s_aux.sof * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.dob_stem[:,irt] .= sen_geo.s_aux.dob * SPECTRA.ρ_STEM;
        sen_geo.auxil.dof_stem[:,irt] .= sen_geo.s_aux.dof * SPECTRA.ρ_STEM;
        sen_geo.auxil.so_stem[:,irt]  .= sen_geo.s_aux.sob * SPECTRA.ρ_STEM;
    end;

    return nothing
end;

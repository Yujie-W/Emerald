# This file contains function to update the leaf pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-26: add function to update the leaf pressure profile
#     2023-Sep-28: compute leaf critical flow here (to use in stomatal optimization model)
#     2023-Oct-16: make e_crit of leaf to be per leaf area (to use with stomatal optimization models)
#     2024-Feb-28: add LAI <= 0 control
#     2024-Jul-24: use spac cache
#
#######################################################################################################################################################################################################
"""

    leaf_pressure_profile!(config::SPACConfiguration{FT}, leaf::Leaf{FT}, cache::SPACCache{FT}, p_dos::FT) where {FT}

Update the leaf pressure profile, given
- `config` `SPACConfiguration` type struct
- `leaf` `Leaf` type struct
- `cache` `SPACCache` type struct
- `p_dos` pressure at the dosing point `[MPa]`

"""
function leaf_pressure_profile!(config::SPACConfiguration{FT}, leaf::Leaf{FT}, cache::SPACCache{FT}, p_dos::FT) where {FT}
    if leaf.xylem.trait.area <= 0
        return nothing
    end;

    # run the pressure profile calculation only if xylem area > 0
    leaf.xylem.auxil.pressure[1] = p_dos;
    leaf.xylem.auxil.e_crit = critical_flow(config, leaf.xylem, cache, leaf.energy.s_aux.t, leaf.xylem.auxil.e_crit) / leaf.xylem.trait.area;
    xylem_pressure_profile!(leaf.xylem, leaf.energy.s_aux.t);
    extraxylary_pressure_profile!(leaf);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function leaf_pressure_profiles!
#     2024-Feb-28: add LAI <= 0 control
#
#######################################################################################################################################################################################################
"""

    leaf_pressure_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Set up leaf pressure profile for each leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `BulkSPAC` type struct

"""
function leaf_pressure_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    # run the pressure profile calculation for each leaf layer only if LAI > 0
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;

    for i in eachindex(branches)
        leaf_pressure_profile!(config, leaves[i], spac.cache, branches[i].xylem.auxil.pressure[end]);
    end;

    return nothing
end;

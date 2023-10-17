# This file contains function to update the leaf pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-26: add function to update the leaf pressure profile
#     2023-Sep-28: compute leaf critical flow here (to use in stomatal optimization model)
#     2023-Oct-16: make e_crit of leaf to be per leaf area (to use with stomatal optimization models)
#
#######################################################################################################################################################################################################
"""

    leaf_pressure_profile!(config::SPACConfiguration{FT}, leaf::Leaf{FT}, p_dos::FT) where {FT}

Update the leaf pressure profile, given
- `config` `SPACConfiguration` type struct
- `leaf` `Leaf` type struct
- `p_dos` pressure at the dosing point `[MPa]`

"""
function leaf_pressure_profile!(config::SPACConfiguration{FT}, leaf::Leaf{FT}, p_dos::FT) where {FT}
    leaf.xylem.auxil.pressure[1] = p_dos;
    xylem_pressure_profile!(leaf.xylem, leaf.energy.auxil.t);
    extraxylary_pressure_profile!(leaf);

    leaf.xylem.auxil.e_crit = critical_flow(config, leaf.xylem, leaf.energy.auxil.t, leaf.xylem.auxil.e_crit) / leaf.xylem.state.area;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function leaf_pressure_profiles!
#
#######################################################################################################################################################################################################
"""

    leaf_pressure_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Set up leaf pressure profile for each leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `BulkSPAC` type struct

"""
function leaf_pressure_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;

    for i in eachindex(branches)
        leaf_pressure_profile!(config, leaves[i], branches[i].xylem.auxil.pressure[end]);
    end;

    return nothing
end;

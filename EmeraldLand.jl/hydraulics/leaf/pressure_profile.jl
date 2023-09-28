# This file contains function to update the leaf pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-26: add function to update the leaf pressure profile
#
#######################################################################################################################################################################################################
"""

    leaf_pressure_profile!(leaf::Leaf{FT}, p_dos::FT) where {FT}

Update the leaf pressure profile, given
- `leaf` `Leaf` type struct
- `p_dos` pressure at the dosing point `[MPa]`

"""
function leaf_pressure_profile!(leaf::Leaf{FT}, p_dos::FT) where {FT}
    leaf.xylem.auxil.pressure[1] = p_dos;
    xylem_pressure_profile!(leaf.xylem, leaf.energy.auxil.t);
    extraxylary_pressure_profile!(leaf);

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

    leaf_pressure_profiles!(spac::MultiLayerSPAC{FT}) where {FT}

Set up leaf pressure profile for each leaf, given
- `spac` `MultiLayerSPAC` type struct

"""
function leaf_pressure_profiles!(spac::MultiLayerSPAC{FT}) where {FT}
    (; BRANCHES, LEAVES) = spac;

    for i in eachindex(BRANCHES)
        leaf_pressure_profile!(LEAVES[i], BRANCHES[i].xylem.auxil.pressure[end]);
    end;

    return nothing
end;

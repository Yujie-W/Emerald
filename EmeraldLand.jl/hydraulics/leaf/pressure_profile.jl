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

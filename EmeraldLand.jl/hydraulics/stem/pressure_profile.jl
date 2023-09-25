# This file contains functions related to stem pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-25: add function stem_pressure_profile!
#
#######################################################################################################################################################################################################
"""

    stem_pressure_profile!(stem::Stem{FT}, p_dos::FT) where {FT}

Update the stem xylem pressure profile, given
- `stem` `Stem` type struct
- `p_dos` Downstream xylem pressure

"""
function stem_pressure_profile!(stem::Stem{FT}, p_dos::FT) where {FT}
    stem.xylem.auxil.pressure[1] = p_dos;
    xylem_pressure_profile!(stem.xylem, stem.energy.auxil.t);

    return nothing
end;

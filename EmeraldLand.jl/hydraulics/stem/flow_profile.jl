# This file contains functions related to stem flow profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-5: add function stem_flow_profile!
#
#######################################################################################################################################################################################################
"""

    stem_flow_profile!(stem::Stem{FT}, flow::FT) where {FT}

Set up stem flow profile, given
- `stem` `Stem` type struct
- `flow` Flow rate out of the stem

"""
function stem_flow_profile!(stem::Stem{FT}, flow::FT) where {FT}
    set_flow_profile!(stem.xylem, flow);

    return nothing
end

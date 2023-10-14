# This file contains function to read beta factor out from the stomatal models

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-11: rename the methods from β_factor to read_β
#
#######################################################################################################################################################################################################
"""

    read_β(leaf::Leaf{FT}) where {FT}

Return the β factor, given
- `leaf` Leaf type

"""
function read_β(leaf::Leaf{FT}) where {FT}
    return leaf.flux.auxil.β
end;

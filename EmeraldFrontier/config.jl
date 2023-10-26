#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Apr-13: add function to create spac configuration
#
#######################################################################################################################################################################################################
"""

    spac_config(gmdict::Dict)

Create a SPAC configuration struct, given
- `gmdict` Dictionary of GriddingMachine data in a grid

"""
function spac_config(gmdict::Dict)
    FT = gmdict["FT"];

    return SPACConfiguration{FT}(MESSAGE_LEVEL = gmdict["MESSAGE_LEVEL"])
end;

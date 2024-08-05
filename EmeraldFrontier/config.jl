#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Apr-13: add function to create spac configuration
#
#######################################################################################################################################################################################################
"""

    spac_config(gm_dict::Dict)

Create a SPAC configuration struct, given
- `gm_dict` Dictionary of GriddingMachine data in a grid

"""
function spac_config(gm_dict::Dict)
    config = SPACConfiguration(gm_dict["FT"]);
    config.MESSAGE_LEVEL = gm_dict["MESSAGE_LEVEL"];
    config.DIM_PPAR_BINS = 10;

    return config
end;

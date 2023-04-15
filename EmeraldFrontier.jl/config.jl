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
    FT = gm_dict["FT"];
    _gcf = GeneralConfiguration();
    _dims = SPACDimension(_gcf, 10);
    _config = SPACConfiguration{FT,_dims}(_gcf);

    return _config
end

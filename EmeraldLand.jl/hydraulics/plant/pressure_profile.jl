# This file contains functions to update the pressure profile along the rhizosphere-root-stem-leaf continuum
# Note that the pressure profile is updated in the following order:
#     1. update the pressure profile of roots
#     2. update the pressure profile of trunk stem
#     3. update the pressure profile of branch stem
#     4. update the pressure profile of leaves


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function plant_pressure_profile!
#
#######################################################################################################################################################################################################
"""

    plant_pressure_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Set up the pressure profile of the plant, given
- `config` `SPACConfiguration` type struct
- `spac` `MultiLayerSPAC` type struct

"""
function plant_pressure_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    root_pressure_profiles!(spac);
    stem_pressure_profiles!(spac);
    leaf_pressure_profiles!(config, spac);

    return nothing
end;

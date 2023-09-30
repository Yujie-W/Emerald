# This file contains functions to update the flow profile along the rhizosphere-root-stem-leaf continuum
# Note that the flow profile is updated in the following order:
#     0. update the flow rates from the buffer system (do this after the updating of the water storage)
#     1. compute the flow rate exiting the leaf (from stomtal conductance)
#     2. compute the flow rate exiting the branch stem (into the leaf)
#     3. compute the flow rate exiting the trunk stem (into the branch stem)
#     4. compute the flow rate exiting the roots (into the root-trunk junction)
# Run this function only after the updating of the water storage

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function plant_flow_profile!
#
#######################################################################################################################################################################################################
"""

    plant_flow_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Set up the flow profile along the rhizosphere-root-stem-leaf continuum, given
- `config` `SPACConfiguration` type struct
- `spac` `MultiLayerSPAC` type struct

"""
function plant_flow_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    leaf_flow_profiles!(config, spac);
    stem_flow_profiles!(spac);
    root_flow_profiles!(config, spac);

    return nothing
end;

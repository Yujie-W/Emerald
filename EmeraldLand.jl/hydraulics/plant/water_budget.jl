# This file contains function to update the water budget of the plant
# Note that the water budget is updated in the following order:
#     1. update the water budget of the leaf water storage
#     2. update the water budget of the branch stem water storage
#     3. update the water budget of the trunk stem water storage
#     4. update the water budget of the root water storage
#     5. update the water budget of the root-trunk junction water storage

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function plant_water_budget!
#
#######################################################################################################################################################################################################
"""

    plant_water_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Set up the water budget of the plant, given
- `spac` `MultiLayerSPAC` type struct
- `δt` time step

"""
function plant_water_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    leaf_water_budgets!(spac, δt);
    stem_water_budgets!(spac, δt);
    root_water_budgets!(spac, δt);
    junction_water_budget!(spac);

    return nothing
end;

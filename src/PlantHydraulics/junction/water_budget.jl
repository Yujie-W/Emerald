# This file contains functions to update the water budget of the root-trunk junction

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function junction_water_budget!
#     2024-Jul-24: use computed ∂w∂t to update the water budget
#
#######################################################################################################################################################################################################
"""

    junction_water_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}

Update the water budget of the root-trunk junction, given
- `spac` `BulkSPAC` type struct
- `δt` Time step

"""
function junction_water_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    spac.plant.junction.state.v_storage += spac.plant.junction.auxil.∂w∂t * δt;

    return nothing
end;

# This file contains function to update the water budget of the root

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#
#######################################################################################################################################################################################################
"""

    root_water_budgets!(spac::BulkSPAC{FT}, δt::FT) where {FT}

Set the flow profile of each root, given
- `spac` `BulkSPAC` type struct
- `δt` time step

"""
function root_water_budgets!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    for root in spac.plant.roots
        xylem_water_budget!(root.xylem, root.xylem.auxil, δt);
    end;

    return nothing
end;

# This file contains function to update the water budget of the root

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#
#######################################################################################################################################################################################################
"""

    root_water_budgets!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Set the flow profile of each root, given
- `spac` `MultiLayerSPAC` type struct
- `δt` time step

"""
function root_water_budgets!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    for root in spac.ROOTS
        xylem_water_budget!((root).NS.xylem, (root).NS.xylem.auxil, (root).NS.energy.auxil.t, δt);
    end;

    return nothing
end;

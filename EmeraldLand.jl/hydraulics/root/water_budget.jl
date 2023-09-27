# This file contains function to update the water budget of the root

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#
#######################################################################################################################################################################################################
"""

    root_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}

Set the flow profile of each root, given
- `spac` `MultiLayerSPAC` type struct
- `Δt` time step

"""
function root_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}
    for root in spac.ROOTS
        xylem_water_budget!(root, (root).NS, (root).NS.xylem.auxil, (root).NS.energy.t, Δt);
    end;

    return nothing
end;

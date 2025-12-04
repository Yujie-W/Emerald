# This file contains function to update the water budget of the trunk and branches

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budgets! function
#     2024-Feb-28: add LAI <= 0 control
#
#######################################################################################################################################################################################################
"""

    stem_water_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}

Set the flow profile of each stem (trunk and branches), given
- `spac` `BulkSPAC` type struct
- `δt` time step

"""
function stem_water_budgets!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    # run the water budget for each stem layer only if LAI > 0
    xylem_water_budget!(spac.plant.trunk.xylem, spac.plant.trunk.xylem.auxil, δt);

    for stem in spac.plant.branches
        xylem_water_budget!(stem.xylem, stem.xylem.auxil, δt);
    end;

    return nothing
end;

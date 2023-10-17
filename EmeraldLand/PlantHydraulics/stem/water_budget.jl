# This file contains function to update the water budget of the trunk and branches

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budgets! function
#
#######################################################################################################################################################################################################
"""

    stem_water_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}

Set the flow profile of each stem (trunk and branches), given
- `spac` `BulkSPAC` type struct
- `δt` time step

"""
function stem_water_budgets!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    xylem_water_budget!(spac.TRUNK.xylem, spac.TRUNK.xylem.auxil, δt);

    for stem in spac.BRANCHES
        xylem_water_budget!(stem.xylem, stem.xylem.auxil, δt);
    end;

    return nothing
end;

# This file contains function to update the water budget of the trunk and branches

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budgets! function
#
#######################################################################################################################################################################################################
"""

    stem_water_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Set the flow profile of each stem (trunk and branches), given
- `spac` `MultiLayerSPAC` type struct
- `δt` time step

"""
function stem_water_budgets!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    xylem_water_budget!(spac.TRUNK, (spac.TRUNK).NS.xylem, (spac.TRUNK).NS.xylem.auxil, (spac.TRUNK).NS.energy.auxil.t, δt);

    for stem in spac.BRANCHES
        xylem_water_budget!(stem, (stem).NS.xylem, (stem).NS.xylem.auxil, (stem).NS.energy.auxil.t, δt);
    end;

    return nothing
end;

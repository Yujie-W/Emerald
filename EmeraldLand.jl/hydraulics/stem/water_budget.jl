# This file contains function to update the water budget of the trunk and branches

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#
#######################################################################################################################################################################################################
"""

    stem_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}

Set the flow profile of each stem (trunk and branches), given
- `spac` `MultiLayerSPAC` type struct
- `Δt` time step

"""
function stem_water_budget(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}
    xylem_water_budget!(spac.TRUNK, (spac.TRUNK).NS, (spac.TRUNK).NS.xylem.auxil, (spac.TRUNK).NS.energy.t, Δt);

    for stem in spac.BRANCHES
        xylem_water_budget!(stem, (stem).NS, (stem).NS.xylem.auxil, (stem).NS.energy.t, Δt);
    end;

    return nothing
end;

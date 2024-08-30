# This file determines what happens when a plant dies

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-30: add function to kill the plant
#
#######################################################################################################################################################################################################
"""

    plant_death!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Kill the plant if the carbon pool is empty, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC

"""
function plant_death!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # if the carbon pool runs out, the plant dies
    if spac.plant.pool.c_pool > 0
        return nothing
    end

    # if the plant is dead, set everything to zero
    # TODO: mass conservation (including the carbon pool)
    shed_leaves!(config, spac);

    # set all xylem area to zero
    for r in spac.plant.roots
        r.xylem.trait.area = 0;
        r.xylem.state.asap = 0;
    end;
    spac.plant.trunk.xylem.trait.area = 0;
    spac.plant.trunk.xylem.state.asap = 0;
    for s in spac.plant.branches
        s.xylem.trait.area = 0;
        s.xylem.state.asap = 0;
    end;

    return nothing
end;

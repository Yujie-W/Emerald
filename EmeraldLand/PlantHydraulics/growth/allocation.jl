# This file contains the algorithm to allocate carbon to different plant organs

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-30: add function to allocate carbon to different plant organs
#
#######################################################################################################################################################################################################
"""

    plant_growth!(spac::BulkSPAC{FT})

Allocate carbon to different plant organs, given
- `spac` `BulkSPAC` structure

"""
function plant_growth!(spac::BulkSPAC{FT}) where {FT}
    plant = spac.plant;

    # determine how much carbon is available for xylem growth
    #     - if not regrow leaf mode, use the max pool threshold
    #     - if regrow leaf mode, use the min pool threshold (immediately after leaf regrow)
    if plant._leaf_regrow
        c_mol = plant.pool.c_pool - plant.pool.c_pool_min;
    else
        c_mol = plant.pool.c_pool - plant.pool_c_pool_max;
    end;

    # if c_mol > 0, allocate the carbon to the xylem
    if c_mol <= 0
        return nothing;
    end;

    # compute the new trunk area to grow
    branches = plant.branches;
    trunk = plant.trunk;
    roots = plant.roots;
    denom = 0;
    for s in branches
        denom += s.xylem.trait.area / (trunk).xylem.trait.area * s.xylem.trait.l * s.xylem.trait.ρ;
    end;
    denom += (trunk).xylem.trait.l * (trunk).xylem.trait.ρ;
    for r in roots
        denom += r.xylem.trait.area / (trunk).xylem.trait.area * r.xylem.trait.l * r.xylem.trait.ρ;
    end;
    denom *= 1000 / 30;
    delta_a = c_mol / denom;

    # compute the energy cost for the new growth to each organ
    for s in branches
        c_mol_s = delta_a * (s).xylem.trait.area / (trunk).xylem.trait.area * (s).xylem.trait.l * (s).xylem.trait.ρ * 1000 / 30;
        allocate_carbon!(s, c_mol_s);
    end;
    c_mol_t = delta_a * (trunk).xylem.trait.l * (trunk).xylem.trait.ρ * 1000 / 30;
    allocate_carbon!(trunk, c_mol_t);
    for r in roots
        c_mol_r = delta_a * (r).xylem.trait.area / (trunk).xylem.trait.area * (r).xylem.trait.l * (r).xylem.trait.ρ * 1000 / 30;
        allocate_carbon!(r, c_mol_r);
    end;

    # update the carbon pool
    plant.pool.c_pool -= c_mol;

    return nothing
end;

allocate_carbon!(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT) where {FT} = recovery_or_growth(organ, c_mol) ? xylem_recovery!(organ, c_mol) : xylem_growth!(organ, c_mol);

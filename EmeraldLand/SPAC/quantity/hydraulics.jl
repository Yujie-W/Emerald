#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-05: add function to return the plant hydraulic conductance of the whole plant
#
#######################################################################################################################################################################################################
"""

    K_PLANT(spac::BulkSPAC{FT}; include_leaf::Bool = true) where {FT}

Return the plant hydraulic conductance of the whole plant, given
- `spac` `BulkSPAC` SPAC
- `include_leaf` `Bool` include leaf conductance or not

"""
function K_PLANT end;

K_PLANT(spac::BulkSPAC{FT}; include_leaf::Bool = true) where {FT} = (
    plant = spac.plant;

    r_plant::FT = 0;
    # iterate through the roots to get the root k
    k_root::FT = 0;
    for r in plant.roots
        k_root += xylem_conductance(r);
    end;
    r_plant += 1 / k_root;

    r_plant += 1 / xylem_conductance(plant.trunk);

    # iterate through the stem and leaf to get the stem and leaf k
    k_canopy::FT = 0;
    for i in eachindex(plant.leaves)
        r_stem = 1 / xylem_conductance(plant.branches[i]);
        r_leaf = include_leaf ? 1 / xylem_conductance(plant.leaves[i]) : FT(0);
        k_canopy += 1 / (r_stem + r_leaf);
    end;
    r_plant += 1 / k_canopy;

    return 1 / r_plant
);

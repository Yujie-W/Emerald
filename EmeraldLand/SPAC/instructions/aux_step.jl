# This file contains functions to update the auxiliary variables at big time step

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_step_auxils! function
#     2024-Jan-24: update leaf boundary layer conductance based on wind speed and leaf width
#     2024-Jul-24: add leaf shedded flag
#
#######################################################################################################################################################################################################
"""

    step_aux!(spac::BulkSPAC{FT}) where {FT}

Update the auxiliary variables at big time step, given
- `spac` `BulkSPAC` SPAC

"""
function step_aux! end;

step_aux!(spac::BulkSPAC{FT}) where {FT} = (
    if spac.plant._leaf_shedded
        return nothing
    end;

    airs = spac.airs;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;

    # clear the integrated values
    for leaf in leaves
        step_aux!(leaf);
    end;

    # update leaf boundary layer conductance
    for i in eachindex(leaves)
        leaf = leaves[i];
        air = airs[lindex[i]];
        leaf.flux.auxil.g_CO₂_b = FT(0.14) * sqrt(air.auxil.wind / (FT(0.72) * leaf.bio.trait.width));
    end;

    return nothing
);

step_aux!(leaf::Leaf{FT}) where {FT} = (
    leaf.flux.auxil.∫∂w∂t_out = 0;

    return nothing
);

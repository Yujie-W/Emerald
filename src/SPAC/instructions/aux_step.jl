# This file contains functions to update the auxiliary variables at big time step

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_step_auxils! function
#     2024-Jan-24: update leaf boundary layer conductance based on wind speed and leaf width
#     2024-Jul-24: add leaf shedded flag
#     2024-Jul-30: compute OCS boundary layer conductance as well
#     2024-Nov-05: remove leaf shedded flag
#
#######################################################################################################################################################################################################
"""

    step_aux!(spac::BulkSPAC{FT}) where {FT}

Update the auxiliary variables at big time step, given
- `spac` `BulkSPAC` SPAC

"""
function step_aux! end;

step_aux!(spac::BulkSPAC{FT}) where {FT} = (
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
        leaf.flux.auxil.g_OCS_b = leaf.flux.auxil.g_CO₂_b / FT(1.21);
    end;

    return nothing
);

step_aux!(leaf::CanopyLayer{FT}) where {FT} = (
    leaf.flux.auxil.∫∂c∂t_in = 0;
    leaf.flux.auxil.∫∂w∂t_out = 0;

    return nothing
);

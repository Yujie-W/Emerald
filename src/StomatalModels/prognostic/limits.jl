# This file contains functions to limit stomatal conductance within its structural limits

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from gsw_control! to limit_stomatal_conductance!
#     2024-Oct-30: add leaf connection check
#
#######################################################################################################################################################################################################
"""

    limit_stomatal_conductance!(leaf::Union{CanopyLayer{FT}, Leaf{FT}}) where {FT}

Limit stomatal conductance for H₂O for
- `leaf` `Leaf` type struct

"""
function limit_stomatal_conductance! end;

limit_stomatal_conductance!(leaf::CanopyLayer{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.state.g_H₂O_s .= 0;

        return nothing
    end;

    # update the g_H₂O_s based on the structural limits if leaf xylem is connected
    f_dif = relative_diffusive_coefficient(leaf.energy.s_aux.t);
    g_min = leaf.flux.trait.g_limits[1] * f_dif;
    g_max = leaf.flux.trait.g_limits[2] * f_dif;

    @. leaf.flux.state.g_H₂O_s = max(g_min, min(g_max, leaf.flux.state.g_H₂O_s));

    return nothing
);

limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.state.g_H₂O_s = 0;

        return nothing
    end;

    # update the g_H₂O_s based on the structural limits if leaf xylem is connected
    f_dif = relative_diffusive_coefficient(leaf.energy.s_aux.t);
    g_min = leaf.flux.trait.g_limits[1] * f_dif;
    g_max = leaf.flux.trait.g_limits[2] * f_dif;

    leaf.flux.state.g_H₂O_s = max(g_min, min(g_max, leaf.flux.state.g_H₂O_s));

    return nothing
);

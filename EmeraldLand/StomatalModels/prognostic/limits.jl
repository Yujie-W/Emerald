# This file contains functions to limit stomatal conductance within its structural limits

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from gsw_control! to limit_stomatal_conductance!
#
#######################################################################################################################################################################################################
"""

    limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT}

Limit stomatal conductance for H₂O for
- `leaf` `Leaf` type struct

"""
function limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT}
    f_dif = relative_diffusive_coefficient(leaf.energy.s_aux.t);
    g_min = leaf.flux.state.g_limits[1] * f_dif;
    g_max = leaf.flux.state.g_limits[2] * f_dif;

    # for sunlit leaf
    for i in eachindex(leaf.flux.state.g_H₂O_s_sunlit)
        if leaf.flux.state.g_H₂O_s_sunlit[i] < g_min
            leaf.flux.state.g_H₂O_s_sunlit[i] = g_min;
        end;
        if leaf.flux.state.g_H₂O_s_sunlit[i] > g_max
            leaf.flux.state.g_H₂O_s_sunlit[i] = g_max;
        end;
    end;

    # for shaded leaf
    if leaf.flux.state.g_H₂O_s_shaded < g_min
        leaf.flux.state.g_H₂O_s_shaded = g_min;
    end;
    if leaf.flux.state.g_H₂O_s_shaded > g_max
        leaf.flux.state.g_H₂O_s_shaded = g_max;
    end;

    return nothing
end;

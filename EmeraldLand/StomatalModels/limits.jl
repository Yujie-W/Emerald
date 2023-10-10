#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from gsw_control! to limit_stomatal_conductance!
#     2022-Jul-01: add method for Leaf
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""

    limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT}

Limit stomatal conductance for H₂O for
- `leaf` `Leaf` type struct

"""
function limit_stomatal_conductance! end;

limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT} = (
    _ratio = relative_diffusive_coefficient(leaf.energy.auxil.t);
    _g_min = leaf.flux.state.g_limits[1] * _ratio;
    _g_max = leaf.flux.state.g_limits[2] * _ratio;

    # for sunlit leaf
    for i in eachindex(leaf.flux.state.g_H₂O_s_sunlit)
        if leaf.flux.state.g_H₂O_s_sunlit[i] < _g_min
            leaf.flux.state.g_H₂O_s_sunlit[i] = _g_min
        end;
        if leaf.flux.state.g_H₂O_s_sunlit[i] > _g_max
            leaf.flux.state.g_H₂O_s_sunlit[i] = _g_max
        end;
    end;

    # for shaded leaf
    if leaf.flux.state.g_H₂O_s_shaded < _g_min
        leaf.flux.state.g_H₂O_s_shaded = _g_min
    end;
    if leaf.flux.state.g_H₂O_s_shaded > _g_max
        leaf.flux.state.g_H₂O_s_shaded = _g_max
    end;

    return nothing
);

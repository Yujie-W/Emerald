#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from gsw_control! to limit_stomatal_conductance!
#     2022-Jul-01: add method for Leaf
#     2022-Jul-01: add method for Leaves2D
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""

    limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT}
    limit_stomatal_conductance!(leaves::Leaves2D{FT}) where {FT}

Limit stomatal conductance for H₂O for
- `leaf` `Leaf` type struct
- `leaves` `Leaves2D` type struct

"""
function limit_stomatal_conductance! end

limit_stomatal_conductance!(leaves::Leaves2D{FT}) where {FT} = (
    (; G_LIMITS) = leaves;

    _ratio = relative_diffusive_coefficient(leaves.t);
    _g_min = G_LIMITS[1] * _ratio;
    _g_max = G_LIMITS[2] * _ratio;

    # for sunlit leaves
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        if leaves.g_H₂O_s_sunlit[_i] < _g_min
            leaves.g_H₂O_s_sunlit[_i] = _g_min
        end;
        if leaves.g_H₂O_s_sunlit[_i] > _g_max
            leaves.g_H₂O_s_sunlit[_i] = _g_max
        end;
    end;

    # for shaded leaves
    if leaves.g_H₂O_s_shaded < _g_min
        leaves.g_H₂O_s_shaded = _g_min
    end;
    if leaves.g_H₂O_s_shaded > _g_max
        leaves.g_H₂O_s_shaded = _g_max
    end;

    return nothing
);

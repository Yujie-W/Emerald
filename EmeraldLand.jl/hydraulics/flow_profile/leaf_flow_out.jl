#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: add function from xylem_flow_profile! to set_leaf_flow! to be more specific in function name
#
#######################################################################################################################################################################################################
"""

    set_leaf_flow_out!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat}
    set_leaf_flow_out!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Set the flow out from each leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type struct

"""
function set_leaf_flow_out! end

set_leaf_flow_out!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat} = (
    (; ALLOW_LEAF_CONDENSATION) = config;
    (; AIR, LEAF) = spac;

    _g = 1 / (1 / LEAF.g_H₂O_s + 1 / (FT(1.35) * LEAF.g_CO₂_b));
    _d = saturation_vapor_pressure(LEAF.t, LEAF.HS.p_leaf * 1000000) - AIR.p_H₂O;
    ALLOW_LEAF_CONDENSATION ? nothing : _d = max(_d, 0);
    _f = _g * _d / AIR.P_AIR;
    set_flow_out!(LEAF.HS.FLOW, _f);

    return nothing
);

set_leaf_flow_out!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; ALLOW_LEAF_CONDENSATION, DIM_LAYER) = config;
    (; AIR, CANOPY, LEAVES, LEAVES_INDEX) = spac;

    for _i in eachindex(LEAVES)
        _p_sl = CANOPY.OPTICS.p_sunlit[DIM_LAYER + 1 - _i];

        _g_sh = 1 / (1 /LEAVES[_i].g_H₂O_s_shaded + 1 / (FT(1.35) * LEAVES[_i].g_CO₂_b));
        _g_sl = 0;
        for _j in eachindex(LEAVES[_i].g_H₂O_s_sunlit)
            _g_sl += 1 / (1 /LEAVES[_i].g_H₂O_s_sunlit[_j] + 1 / (FT(1.35) * LEAVES[_i].g_CO₂_b));
        end;
        _g_sl /= length(LEAVES[_i].g_H₂O_s_sunlit);

        _g = _g_sh * (1 - _p_sl) + _g_sl * _p_sl;
        _d = saturation_vapor_pressure(LEAVES[_i].t, LEAVES[_i].HS.p_leaf * 1000000) - AIR[LEAVES_INDEX[_i]].p_H₂O;
        ALLOW_LEAF_CONDENSATION ? nothing : _d = max(_d, 0);
        _f = _g * _d / AIR[LEAVES_INDEX[_i]].P_AIR;

        set_flow_out!(LEAVES[_i].HS.FLOW, _f);
    end;

    return nothing
);

# This file contains functions related to leaf flow profile

flow_in(leaf::Leaf{FT}) where {FT} = flow_in(leaf.xylem);

flow_in(leaf::Leaves2D{FT}) where {FT} = flow_in(leaf.NS);

flow_out(leaf::Leaf{FT}) where {FT} = flow_out(leaf.xylem) + leaf.capacitor.auxil.flow;

flow_out(leaf::Leaves2D{FT}) where {FT} = flow_out(leaf.NS);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: add function leaf_flow_profiles!
#     2023-Sep-28: account for the buffer from capacitor when running under non-steady state mode
#
#######################################################################################################################################################################################################
"""

    leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Set the flow out from each leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `MultiLayerSPAC` type struct

"""
function leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    # compute the flow rate exiting the leaf based on sunlit and shaded fractions and update it to the leaf of a MultiLayerSPAC
    #     LEAVES index is from lower to upper, and thus the sunlit leaves fraction is DIM_LAYER + 1 - i
    #     AIR index is also from lower to upper, but there are some layers are used by trunk so that it need to be indexed through LEAVES_INDEX

    (; ALLOW_LEAF_CONDENSATION, DIM_LAYER) = config;
    (; AIR, CANOPY, LEAVES, LEAVES_INDEX) = spac;

    for i in eachindex(LEAVES)
        f_sl = CANOPY.OPTICS.p_sunlit[DIM_LAYER + 1 - i];

        g_sh = 1 / (1 /LEAVES[i].g_H₂O_s_shaded + 1 / (FT(1.35) * LEAVES[i].g_CO₂_b));
        g_sl = 0;
        for j in eachindex(LEAVES[i].g_H₂O_s_sunlit)
            g_sl += 1 / (1 / LEAVES[i].g_H₂O_s_sunlit[j] + 1 / (FT(1.35) * LEAVES[i].g_CO₂_b));
        end;
        g_sl /= length(LEAVES[i].g_H₂O_s_sunlit);

        g = g_sh * (1 - f_sl) + g_sl * f_sl;
        d = saturation_vapor_pressure(LEAVES[i].NS.energy.auxil.t, LEAVES[i].NS.capacitor.auxil.p_leaf * 1000000) - AIR[LEAVES_INDEX[i]].p_H₂O;
        ALLOW_LEAF_CONDENSATION ? nothing : d = max(d, 0);
        f = g * d / AIR[LEAVES_INDEX[i]].P_AIR;

        # set_flow_out!(LEAVES[i].HS.FLOW, f);
        set_flow_profile!(LEAVES[i].NS.xylem, f - LEAVES[i].NS.capacitor.auxil.flow);
    end;

    return nothing
end;
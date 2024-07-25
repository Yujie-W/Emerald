# This file contains the functions to compute the ∂Θ∂E for WangSM optimality model

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air;

    return leaf.flux.auxil.a_n_shaded / max(eps(FT), (leaf.xylem.auxil.e_crit - e))
);

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air;

    return leaf.flux.auxil.a_n_sunlit[ind] / max(eps(FT), (leaf.xylem.auxil.e_crit - e))
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-24: add function to update the ∂Θ∂E for sunlit leaves (matrix version)
#     2024-Jul-25: compute ∂Θ∂E for the entire leaf layer (no sunlit and shaded fractions as the sunlit and shaded part could be on the same leaf)
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E!(sm::WangSM{FT}, leaf::Leaf{FT}) where {FT}

Update the ∂Θ∂E for sunlit leaves, given
- `sm` `WangSM` type WangSM
- `leaf` `Leaf` type leaf

"""
∂Θ∂E!(sm::WangSM{FT}, leaf::Leaf{FT}) where {FT} = (
    leaf.flux.auxil.∂Θ∂E = leaf.flux.auxil.a_n / max(eps(FT), (leaf.xylem.auxil.e_crit - flow_out(leaf)));

    return nothing
);

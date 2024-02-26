# This file contains the functions to compute the ∂Θ∂E for WangSM optimality model

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air;

    return leaf.flux.auxil.a_n_shaded / max(eps(FT), (leaf.xylem.auxil.e_crit - e))
);

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air;

    return leaf.flux.auxil.a_n_sunlit[ind] / max(eps(FT), (leaf.xylem.auxil.e_crit - e))
);

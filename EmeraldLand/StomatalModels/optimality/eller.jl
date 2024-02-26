# This file contains function to compute the partial derivatives of Θ versus E for EllerSM optimality model

∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.state.area;

    dedp1 = ∂E∂P(leaf, e; δe = δe);
    dedp2 = ∂E∂P(leaf, e; δe = -δe);
    dkde  = (dedp2 - dedp1) / δe;

    return dkde * leaf.flux.auxil.a_n_shaded / dedp1
);

∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.state.area;

    dedp1 = ∂E∂P(leaf, e; δe = δe);
    dedp2 = ∂E∂P(leaf, e; δe = -δe);
    dkde  = (dedp2 - dedp1) / δe;

    return dkde * leaf.flux.auxil.a_n_sunlit[ind] / dedp1
);

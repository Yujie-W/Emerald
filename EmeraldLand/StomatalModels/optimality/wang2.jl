# This file contains functions to compute the ∂Θ∂E for the Wang2SM optimization model (a modification of the AndereggSM optimization model)

∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; A) = sm;

    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.trait.area;

    dedp = ∂E∂P(leaf, e; δe = δe) / leaf.xylem.trait.area;

    return (-A * leaf.capacitor.state.p_leaf * leaf.flux.auxil.a_n_shaded) / dedp
);

∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; A) = sm;

    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.trait.area;

    dedp = ∂E∂P(leaf, e; δe = δe) / leaf.xylem.trait.area;

    return (-A * leaf.capacitor.state.p_leaf * leaf.flux.auxil.a_n_sunlit[ind]) / dedp
);

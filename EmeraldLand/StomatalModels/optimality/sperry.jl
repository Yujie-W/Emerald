# This file contains functions to compute the ∂Θ∂E for the SperrySM

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2023-Oct-16: make sure maximum gsc does not exceed g_CO₂_b
#
#######################################################################################################################################################################################################
∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.state.area;

    dEdP_1 = ∂E∂P(leaf, e; δe = δe);
    dEdP_2 = ∂E∂P(leaf, e; δe = -δe);
    dEdP_m = ∂E∂P(leaf, FT(0); δe = δe);
    dKdE   = (dEdP_2 - dEdP_1) / δe;

    # compute maximum A
    ghm = leaf.xylem.auxil.e_crit / d * air.state.p_air;
    if ghm < FT(1.35) * leaf.flux.auxil.g_CO₂_b
        gsm = 1 / (1 / ghm - 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
        gcm = 1 / (FT(1.6) / gsm + 1 / leaf.flux.auxil.g_CO₂_b);
    else
        gsm = FT(Inf);
        gcm = leaf.flux.auxil.g_CO₂_b;
    end;
    photosynthesis_only!(leaf.photosystem, air, gcm, leaf.flux.auxil.ppar_shaded, leaf.energy.auxil.t);
    am = leaf.photosystem.auxil.a_n;

    return dKdE * am / dEdP_m
);

∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.state.area;

    dEdP_1 = ∂E∂P(leaf, e; δe = δe);
    dEdP_2 = ∂E∂P(leaf, e; δe = -δe);
    dEdP_m = ∂E∂P(leaf, FT(0); δe = δe);
    dKdE   = (dEdP_2 - dEdP_1) / δe;

    # compute maximum A
    ghm = leaf.xylem.auxil.e_crit / d * air.state.p_air;
    if ghm < FT(1.35) * leaf.flux.auxil.g_CO₂_b
        gsm = 1 / (1 / ghm - 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
        gcm = 1 / (FT(1.6) / gsm + 1 / leaf.flux.auxil.g_CO₂_b);
    else
        gsm = FT(Inf);
        gcm = leaf.flux.auxil.g_CO₂_b;
    end;
    photosynthesis_only!(leaf.photosystem, air, gcm, leaf.flux.auxil.ppar_sunlit[ind], leaf.energy.auxil.t);
    am = leaf.photosystem.auxil.a_n;

    return dKdE * am / dEdP_m
);

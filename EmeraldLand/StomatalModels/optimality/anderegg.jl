# This file contains functions to compute the partial derivatives of Θ versus E for AndereggSM optimality model

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-11: add method for AndereggSM model
#     2022-Jul-11: add method for EllerSM model
#     2022-Jul-11: add method for SperrySM model
#     2022-Jul-11: add method for WangSM model
#     2022-Jul-11: add method for Wang2SM model
#     2023-Mar-02: add a eps(FT) controller to (e_crit - e)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P
- `ind` Index of sunlit leaf

"""
∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; A, B) = sm;

    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.state.area;

    dedp = ∂E∂P(leaf, e; δe = δe) / leaf.xylem.state.area;

    return (-2 * A * leaf.capacitor.auxil.p_leaf + B) / dedp
);

∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; A, B) = sm;

    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e  = gh * d / air.state.p_air * leaf.xylem.state.area;

    dedp = ∂E∂P(leaf, e; δe = δe) / leaf.xylem.state.area;

    return (-2 * A * leaf.capacitor.auxil.p_leaf + B) / dedp
);

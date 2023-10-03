#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-07: add the method for Leaf
#     2022-Jul-07: add the method for Leaves2D (shaded leaves)
#     2022-Jul-07: add the method for Leaves2D (sunlit leaves)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Aug-23: add nan check for some methods (not all, to do later)
#     2023-Aug-27: make sure D >= 1 when compute dAdE
# To do
#     TODO: figure out why the ratios are 1.35 and 1.6, and make them more accurate
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT}
    ∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int) where {FT}

Return the partial derivative of A per E, given
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `leaves` `Leaves2D` type leaf
- `ind` Index of the leaves

"""
function ∂A∂E end

∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = (
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s_shaded;
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e1  = _gh1 * _d / P_AIR;
    _a1  = leaves.a_net_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b);
    photosynthesis_only!(leaves.NS.photosystem, air, _gc2, leaves.ppar_shaded, leaves.NS.energy.auxil.t);
    _e2 = _gh2 * _d / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);

∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int) where {FT} = (
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s_sunlit[ind];
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e1  = _gh1 * _d / P_AIR;
    _a1  = leaves.a_net_sunlit[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b);
    photosynthesis_only!(leaves.NS.photosystem, air, _gc2, leaves.ppar_sunlit[ind], leaves.NS.energy.auxil.t);
    _e2 = _gh2 * _d / P_AIR;
    _a2 = leaves.PSM.a_net;

    _dade = (_a2 - _a1) / (_e2 - _e1);
    if isnan(_dade)
        @info "Debugging" _e1 _a1 _gs1 _gh1 _e2 _a2 _gs2 _gh2 _gc2;
        error("NaN detected when computing ∂A∂E!");
    end;

    return _dade
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂R∂E
#
#######################################################################################################################################################################################################
"""

    ∂R∂E(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT}

Returns the marginal increase in leaf respiration rate per transpiration rate, given
- `lf` `Leaf`, `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions

"""
function ∂R∂E end

∂R∂E(lf::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = ∂R∂E(lf.SM, lf, air);

∂R∂E(sm::WangSM{FT}, lf::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = ∂R∂T(lf) * ∂T∂E(lf, air, sm.f_view);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂T∂E
#
#######################################################################################################################################################################################################
"""

    ∂T∂E(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}, f_view::FT) where {FT}

Returns the marginal increase in leaf temperature per transpiration rate, given
- `lf` `Leaf`, `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions
- `f_view` Ratio that leaf area is exposed to external sources/sinks (not other leaves, e.g., 2/LAI for canopy on average)

"""
function ∂T∂E end

∂T∂E(lf::Leaves2D{FT}, air::AirLayer{FT}, f_view::FT) where {FT} = ∂T∂E(lf.BIO, lf, air, f_view);

∂T∂E(bio::HyperLeafBio{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, f_view::FT) where {FT} = ∂T∂E(f_view, leaves.t, leaves.WIDTH, air.wind, 1 - bio.auxil.τ_LW);

∂T∂E(f_view::FT, t::FT, width::FT, wind::FT, ϵ::FT) where {FT} = (
    _λ = latent_heat_vapor(t) * M_H₂O(FT);
    _g = FT(1.4) * FT(0.135) * sqrt(wind / (FT(0.72) * width));
    _d = 2 * CP_D_MOL(FT) * _g + 4 * f_view * K_STEFAN(FT) * ϵ * t ^ 3;

    return _λ / _d
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-11: add docs
#
#######################################################################################################################################################################################################
"""
This function returns the marginal risk for stomatal opening. This function supports a variety of optimality models for
- Leaf
- Leaves2D (ind=NA for shaded leaves, ind>1 for sunlit leaves)

"""
function ∂Θ∂E end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for EllerSM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for SperrySM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for WangSM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for Wang2SM model on Leaves2D for shaded leaves
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; A, B) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaves.a_net_shaded / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, leaves.HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaves, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / _d * P_AIR;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaves.g_CO₂_b);
    photosynthesis_only!(leaves.NS.photosystem, air, _gcm, leaves.ppar_shaded, leaves.t);
    _am = leaves.PSM.a_net;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    return leaves.a_net_shaded / max(eps(FT), (leaves.NS.xylem.auxil.e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; A) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaves.a_net_shaded) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for EllerSM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for SperrySM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for WangSM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for Wang2SM model on Leaves2D for sunlit leaves
#     2023-Mar-02: add a eps(FT) controller to (e_crit - e)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; A, B) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaves.a_net_sunlit[ind] / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaves, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / _d * P_AIR;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaves.g_CO₂_b);
    photosynthesis_only!(leaves.NS.photosystem, air, _gcm, leaves.ppar_sunlit[ind], leaves.t);
    _am = leaves.PSM.a_net;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    return leaves.a_net_sunlit[ind] / max(eps(FT), (leaves.NS.xylem.auxil.e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; A) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * _d / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaves.a_net_sunlit[ind]) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function for nocturnal stomatal conductance
#
#######################################################################################################################################################################################################
"""
This function returns the ∂Θₙ∂E for nocturnal stomatal opening. Currently this function only supports WangSM which has been published for the purpose of computing nocturnal stomatal conductance.
    Supports to other optimality models will be added later when I am ready to test those.

"""
function ∂Θₙ∂E end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for WangSM model on Leaf for nocturnal transpiration
#     2022-Jul-11: add method for WangSM model on Leaves2D for nocturnal transpiration
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θₙ∂E(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT}

Return the ∂Θ∂E for nocturnal stomatal opening, given
- `lf` `Leaf`, `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions

"""
∂Θₙ∂E(lf::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = ∂Θₙ∂E(lf.SM, lf, air);

∂Θₙ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = (
    (; F_FITNESS) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    _p_s = saturation_vapor_pressure(leaves.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.p_H₂O);

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc = 1 / (FT(1.6) / _gs + 1 / leaves.g_CO₂_b);
    _e  = _gh * _d / P_AIR;
    photosynthesis_only!(leaves.NS.photosystem, air, _gc, sm.ppar_mem, leaves.t);
    _a  = leaves.PSM.a_net;

    return _a / max(eps(FT), (HS._e_crit - _e)) * F_FITNESS
);

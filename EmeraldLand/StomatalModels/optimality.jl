#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-07: add the method for Leaf
#     2022-Jul-07: add the method for Leaf (shaded leaf)
#     2022-Jul-07: add the method for Leaf (sunlit leaf)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Aug-23: add nan check for some methods (not all, to do later)
#     2023-Aug-27: make sure D >= 1 when compute dAdE
# To do
#     TODO: figure out why the ratios are 1.35 and 1.6, and make them more accurate
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}
    ∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int) where {FT}

Return the partial derivative of A per E, given
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `ind` Index of the leaf

"""
function ∂A∂E end;

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    _p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    _gs1 = leaf.flux.state.g_H₂O_s_shaded;
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e1  = _gh1 * _d / air.state.p_air;
    _a1  = leaf.flux.auxil.a_n_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    photosynthesis_only!(leaf.photosystem, air, _gc2, leaf.flux.auxil.ppar_shaded, leaf.energy.auxil.t);
    _e2 = _gh2 * _d / air.state.p_air;
    _a2 = leaf.photosystem.auxil.a_n;

    return (_a2 - _a1) / (_e2 - _e1)
);

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int) where {FT} = (
    _p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    _gs1 = leaf.flux.state.g_H₂O_s_sunlit[ind];
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e1  = _gh1 * _d / air.state.p_air;
    _a1  = leaf.flux.auxil.a_n_sunlit[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    photosynthesis_only!(leaf.photosystem, air, _gc2, leaf.flux.auxil.ppar_sunlit[ind], leaf.energy.auxil.t);
    _e2 = _gh2 * _d / air.state.p_air;
    _a2 = leaf.photosystem.auxil.a_n;

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

    ∂R∂E(leaf::Union{Leaf{FT}, Leaf{FT}}, air::AirLayer{FT}) where {FT}

Returns the marginal increase in leaf respiration rate per transpiration rate, given
- `leaf` `Leaf`, `Leaf` type leaf
- `air` `AirLayer` type environmental conditions

"""
function ∂R∂E end;

∂R∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = ∂R∂E(leaf.flux.state.stomatal_model, leaf, air);

∂R∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = ∂R∂T(leaf) * ∂T∂E(leaf, air, sm.f_view);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂T∂E
#
#######################################################################################################################################################################################################
"""

    ∂T∂E(leaf::Union{Leaf{FT}, Leaf{FT}}, air::AirLayer{FT}, f_view::FT) where {FT}

Returns the marginal increase in leaf temperature per transpiration rate, given
- `leaf` `Leaf`, `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `f_view` Ratio that leaf area is exposed to external sources/sinks (not other leaf, e.g., 2/LAI for canopy on average)

"""
function ∂T∂E end;

∂T∂E(leaf::Leaf{FT}, air::AirLayer{FT}, f_view::FT) where {FT} = ∂T∂E(leaf.BIO, leaf, air, f_view);

∂T∂E(bio::LeafBio{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, f_view::FT) where {FT} = ∂T∂E(f_view, leaf.t, leaf.bio.width, air.auxil.wind, 1 - bio.auxil.τ_lw);

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
- Leaf (ind=NA for shaded leaf, ind>1 for sunlit leaf)

"""
function ∂Θ∂E end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaf for shaded leaf
#     2022-Jul-11: add method for EllerSM model on Leaf for shaded leaf
#     2022-Jul-11: add method for SperrySM model on Leaf for shaded leaf
#     2022-Jul-11: add method for WangSM model on Leaf for shaded leaf
#     2022-Jul-11: add method for Wang2SM model on Leaf for shaded leaf
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; A, B) = sm;
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P = ∂E∂P(leaf, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P_1 = ∂E∂P(leaf, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaf, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaf.flux.auxil.a_n_shaded / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, leaf.HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P_1 = ∂E∂P(leaf, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaf, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaf, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / _d * air.state.p_air;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaf.flux.auxil.g_CO₂_b);
    photosynthesis_only!(leaf.photosystem, air, _gcm, leaf.flux.auxil.ppar_shaded, leaf.t);
    _am = leaf.photosystem.auxil.a_n;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    _p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    return leaf.flux.auxil.a_n_shaded / max(eps(FT), (leaf.xylem.auxil.e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    (; A) = sm;
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P = ∂E∂P(leaf, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaf.flux.auxil.a_n_shaded) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaf for sunlit leaf
#     2022-Jul-11: add method for EllerSM model on Leaf for sunlit leaf
#     2022-Jul-11: add method for SperrySM model on Leaf for sunlit leaf
#     2022-Jul-11: add method for WangSM model on Leaf for sunlit leaf
#     2022-Jul-11: add method for Wang2SM model on Leaf for sunlit leaf
#     2023-Mar-02: add a eps(FT) controller to (e_crit - e)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}
    ∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; A, B) = sm;
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P = ∂E∂P(leaf, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P_1 = ∂E∂P(leaf, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaf, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaf.flux.auxil.a_n_sunlit[ind] / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P_1 = ∂E∂P(leaf, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaf, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaf, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / _d * air.state.p_air;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaf.flux.auxil.g_CO₂_b);
    photosynthesis_only!(leaf.photosystem, air, _gcm, leaf.flux.auxil.ppar_sunlit[ind], leaf.t);
    _am = leaf.photosystem.auxil.a_n;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    _p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    return leaf.flux.auxil.a_n_sunlit[ind] / max(eps(FT), (leaf.xylem.auxil.e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    (; A) = sm;
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _e  = _gh * _d / air.state.p_air;

    _∂E∂P = ∂E∂P(leaf, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaf.flux.auxil.a_n_sunlit[ind]) / _∂E∂P
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
function ∂Θₙ∂E end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for WangSM model on Leaf for nocturnal transpiration
#     2022-Jul-11: add method for WangSM model on Leaf for nocturnal transpiration
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θₙ∂E(leaf::Union{Leaf{FT}, Leaf{FT}}, air::AirLayer{FT}) where {FT}

Return the ∂Θ∂E for nocturnal stomatal opening, given
- `leaf` `Leaf`, `Leaf` type leaf
- `air` `AirLayer` type environmental conditions

"""
∂Θₙ∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = ∂Θₙ∂E(leaf.flux.state.stomatal_model, leaf, air);

∂Θₙ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    (; F_FITNESS) = sm;
    (; HS) = leaf;

    _p_s = saturation_vapor_pressure(leaf.t, HS.p_leaf * 1000000);
    _d = max(1, _p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    _gs = leaf.flux.state.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    _gc = 1 / (FT(1.6) / _gs + 1 / leaf.flux.auxil.g_CO₂_b);
    _e  = _gh * _d / air.state.p_air;
    photosynthesis_only!(leaf.photosystem, air, _gc, sm.ppar_mem, leaf.t);
    _a  = leaf.photosystem.auxil.a_n;

    return _a / max(eps(FT), (HS._e_crit - _e)) * F_FITNESS
);

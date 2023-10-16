


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

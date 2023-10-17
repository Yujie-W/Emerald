# This file contains functions to compute the ∂Θₙ∂E for nocturnal stomatal opening using WangSM model (no other model supported at the moment)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add method for WangSM model on Leaf for nocturnal transpiration
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    ∂Θₙ∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}

Compute the ∂Θₙ∂E for nocturnal stomatal opening, given
- `leaf` Leaf struct
- `air` AirLayer struct

"""
function ∂Θₙ∂E end;

∂Θₙ∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = ∂Θₙ∂E(leaf.flux.state.stomatal_model, leaf, air);

∂Θₙ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    (; F_FITNESS) = sm;

    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    gs = leaf.flux.state.g_H₂O_s_shaded;
    gh = 1 / (1 / gs + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc = 1 / (FT(1.6) / gs + 1 / leaf.flux.auxil.g_CO₂_b);
    e  = gh * d / air.state.p_air;
    photosynthesis_only!(leaf.photosystem, air, gc, leaf.flux.auxil.ppar_mem, leaf.energy.auxil.t);
    a  = leaf.photosystem.auxil.a_n;

    @show a leaf.xylem.auxil.e_crit e F_FITNESS;

    return a / max(eps(FT), (leaf.xylem.auxil.e_crit - e)) * F_FITNESS
);

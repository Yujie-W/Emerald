# This file contains function to compute the marginal change in nighttime stomatal conductance

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add method for nocturnal transpiration for WangSM model
#     2022-Jul-11: rename function to ∂gₙ∂t
#     2023-Mar-11: limit ∂g∂t within (-0.001, 0.001)
#
#######################################################################################################################################################################################################
"""

    ∂gₙ∂t(leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT}

Return the marginal increase of stomatal conductance, given
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `eff_ϵ` Effective emissivity of the canopy layer (single layer value)

"""
function ∂gₙ∂t end;

∂gₙ∂t(leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = ∂gₙ∂t(leaf.flux.state.stomatal_model, leaf, air, eff_ϵ);

∂gₙ∂t(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = (
    drde = ∂R∂E(leaf, air, eff_ϵ);
    dθde = ∂Θₙ∂E(leaf, air);

    @show drde dθde;

    return max(-0.001, min(0.001, sm.K * (drde - dθde)))
);

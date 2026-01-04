# This file contains function to compute the marginal change in nighttime stomatal conductance

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add method for nocturnal transpiration for WangSM model
#     2023-Mar-11: limit ∂gₙ∂t within (-0.001, 0.001)
#     2023-Oct-25: set the ∂gₙ∂t to -0.001 for other models without nighttime stomatal conductance model
#
#######################################################################################################################################################################################################
"""

    ∂gₙ∂t(leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT}

Return the marginal increase of stomatal conductance, given
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `eff_ϵ` Effective emissivity of the leaf layer (single layer value)

"""
function ∂gₙ∂t end;

∂gₙ∂t(config::SPACConfig{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = ∂gₙ∂t(config.METHODS.STOMATAL_MODEL, leaf, air, eff_ϵ);

∂gₙ∂t(sm::AbstractStomatalConductanceModel{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = FT(-0.001);

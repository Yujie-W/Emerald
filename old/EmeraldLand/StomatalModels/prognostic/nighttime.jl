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

    ∂gₙ∂t(leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}, eff_ϵ::FT) where {FT}

Return the marginal increase of stomatal conductance, given
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `eff_ϵ` Effective emissivity of the canopy layer (single layer value)

"""
function ∂gₙ∂t end;

∂gₙ∂t(leaf::CanopyLayer{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = ∂gₙ∂t(leaf.flux.trait.stomatal_model, leaf, air, eff_ϵ);

∂gₙ∂t(sm::AbstractStomataModel{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = FT(-0.001);

#=
∂gₙ∂t(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = (
    drde = ∂R∂E(leaf, air, eff_ϵ);
    dθde = ∂Θₙ∂E(leaf, air);

    return max(-0.001, min(0.001, sm.K * (drde - dθde)))
);
=#

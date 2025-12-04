# This file contains function to compute the marginal change in daytime stomatal conductance


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#     2022-Jul-11: add method for Leaf shaded leaf
#     2022-Jul-11: add method for Leaf sunlit leaf
#     2023-Mar-11: limit ∂g∂t within (-0.001, 0.001)
#     2023-Oct-17: link beta to that stored in leaf for empirical models
#     2024-Jul-24: add function to update the ∂g∂t for sunlit leaves (matrix version)
#     2024-Jul-25: add support to empirical models
#
#######################################################################################################################################################################################################
"""

    ∂g∂t!(cache::SPACCache{FT}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}) where {FT}

Update the ∂g∂t for sunlit leaves, given
- `cache` `SPACCache` type cache
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type air

"""
function ∂g∂t! end;

∂g∂t!(cache::SPACCache{FT}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = ∂g∂t!(cache, leaf.flux.trait.stomatal_model, leaf, air; δe = δe);

∂g∂t!(cache::SPACCache{FT}, sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaf::CanopyLayer{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    ∂A∂E!(cache, leaf, air);
    ∂Θ∂E!(cache, sm, leaf, air);

    leaf.flux.auxil.∂g∂t .= sm.K .* (leaf.flux.auxil.∂A∂E .- leaf.flux.auxil.∂Θ∂E);

    # set the max and min values
    leaf.flux.auxil.∂g∂t .= max.(-0.001, min.(0.001, leaf.flux.auxil.∂g∂t));

    return nothing
);

∂g∂t!(cache::SPACCache{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::CanopyLayer{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} =
    ∂g∂t!(cache, sm, leaf, air, sm.β.PARAM_Y);

∂g∂t!(cache::SPACCache{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::CanopyLayer{FT}, air::AirLayer{FT}, βt::BetaParameterG1) where {FT} = (
    gsw = empirical_equation(sm, leaf, air; β = leaf.flux.auxil.β);
    leaf.flux.auxil.∂g∂t .= max.(-0.001, min.(0.001, (gsw .- leaf.flux.state.g_H₂O_s) ./ sm.τ));

    return nothing
);

∂g∂t!(cache::SPACCache{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::CanopyLayer{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax) where {FT} = (
    gsw = empirical_equation(sm, leaf, air; β = FT(1));
    leaf.flux.auxil.∂g∂t .= max.(-0.001, min.(0.001, (gsw .- leaf.flux.state.g_H₂O_s) ./ sm.τ));

    return nothing
);

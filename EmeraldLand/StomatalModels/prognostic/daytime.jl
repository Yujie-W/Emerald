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
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT}
    ∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT}

Return the marginal change in stomatal conductance, given
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)
- `ind` Sunlit leaf index within the leaf angular distribution (if present, model is meant for sunlit leaves; otherwise, model is meant for shaded leaves)

"""
function ∂g∂t end;

# for shaded leaves
∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = ∂g∂t(leaf.flux.trait.stomatal_model, leaf, air; δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    return max(-0.001, min(0.001, sm.K * (∂A∂E(leaf, air) - ∂Θ∂E(sm, leaf, air; δe = δe))))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    return ∂g∂t(sm, leaf, air, sm.β.PARAM_Y)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1) where {FT} = (
    gsw = empirical_equation(sm, leaf, air; β = leaf.flux.auxil.β);

    return max(-0.001, min(0.001, (gsw - leaf.flux.state.g_H₂O_s_shaded) / sm.τ))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax) where {FT} = (
    gsw = empirical_equation(sm, leaf, air; β = FT(1));

    return max(-0.001, min(0.001, (gsw - leaf.flux.state.g_H₂O_s_shaded) / sm.τ))
);

# for sunlit leaves
∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = ∂g∂t(leaf.flux.trait.stomatal_model, leaf, air, ind; δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = (
    dade = ∂A∂E(leaf, air, ind);
    dtde = ∂Θ∂E(sm, leaf, air, ind; δe = δe);

    return max(-0.001, min(0.001, sm.K * (dade - dtde)))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT} = ∂g∂t(sm, leaf, air, sm.β.PARAM_Y, ind);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int) where {FT} = (
    gsw = empirical_equation(sm, leaf, air, ind; β = leaf.flux.auxil.β);

    return max(-0.001, min(0.001, (gsw - leaf.flux.state.g_H₂O_s_sunlit[ind]) / sm.τ))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int) where {FT} = (
    gsw = empirical_equation(sm, leaf, air, ind; β = FT(1));

    return max(-0.001, min(0.001, (gsw - leaf.flux.state.g_H₂O_s_sunlit[ind]) / sm.τ))
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-24: add function to update the ∂g∂t for sunlit leaves (matrix version)
#     2024-Jul-25: add support to empirical models
#
#######################################################################################################################################################################################################
"""

    ∂g∂t!(cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}

Update the ∂g∂t for sunlit leaves, given
- `cache` `SPACCache` type cache
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type air

"""
function ∂g∂t! end;

∂g∂t!(cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = ∂g∂t!(cache, leaf.flux.trait.stomatal_model, leaf, air; δe = δe);

∂g∂t!(cache::SPACCache{FT}, sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} = (
    ∂A∂E!(cache, leaf, air);
    ∂Θ∂E!(sm, leaf);

    leaf.flux.auxil.∂g∂t_sunlit .= sm.K .* (leaf.flux.auxil.∂A∂E_sunlit .- leaf.flux.auxil.∂Θ∂E);

    # set the max and min values
    for i in eachindex(leaf.flux.auxil.∂g∂t_sunlit)
        leaf.flux.auxil.∂g∂t_sunlit[i] = max(-0.001, min(0.001, leaf.flux.auxil.∂g∂t_sunlit[i]));
    end;

    return nothing
);

∂g∂t!(cache::SPACCache{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT} =
    ∂g∂t!(cache, sm, leaf, air, sm.β.PARAM_Y);

∂g∂t!(cache::SPACCache{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1) where {FT} = (
    for ind in eachindex(leaf.flux.auxil.∂g∂t_sunlit)
        gsw = empirical_equation(sm, leaf, air, ind; β = leaf.flux.auxil.β);
        leaf.flux.auxil.∂g∂t_sunlit[ind] = max(-0.001, min(0.001, (gsw - leaf.flux.state.g_H₂O_s_sunlit[ind]) / sm.τ));
    end;

    return nothing
);

∂g∂t!(cache::SPACCache{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax) where {FT} = (
    for ind in eachindex(leaf.flux.auxil.∂g∂t_sunlit)
        gsw = empirical_equation(sm, leaf, air, ind; β = FT(1));
        leaf.flux.auxil.∂g∂t_sunlit[ind] = max(-0.001, min(0.001, (gsw - leaf.flux.state.g_H₂O_s_sunlit[ind]) / sm.τ));
    end;

    return nothing
);

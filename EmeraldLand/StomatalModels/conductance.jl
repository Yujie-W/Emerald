#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance change slope. Supported functionalities are
- Leaf
- Leaves2D (ind=NA for shaded, ind>1 for sunlit leaves)

"""
function ∂g∂t end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for Leaves2D shaded leaves
#     2023-Mar-11: limit ∂g∂t within (-0.001, 0.001)
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT}

Return the marginal increase of stomatal conductance, given
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor (only used for empirical models)
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)

"""
∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT} = ∂g∂t(leaves.SM, leaves, air; β = β, δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT} = (
    return max(-0.001, min(0.001, sm.K * (∂A∂E(leaves, air) - ∂Θ∂E(sm, leaves, air; δe = δe))))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT} = (
    return ∂g∂t(sm, leaves, air, sm.β.PARAM_Y; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT} = (
    _gsw = empirical_equation(sm, leaves, air; β = β);

    return max(-0.001, min(0.001, (_gsw - leaves.g_H₂O_s_shaded) / sm.τ))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT} = (
    _gsw = empirical_equation(sm, leaves, air; β = FT(1));

    return max(-0.001, min(0.001, (_gsw - leaves.g_H₂O_s_shaded) / sm.τ))
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for Leaves2D sunlit leaves
#     2023-Mar-11: limit ∂g∂t within (-0.001, 0.001)
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT}

Return the marginal increase of stomatal conductance, given
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Sunlit leaf index within the leaf angular distribution
- `β` Tuning factor (only used for empirical models)
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)

"""
∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT} = ∂g∂t(leaves.SM, leaves, air, ind; β = β, δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT} = (
    _dade = ∂A∂E(leaves, air, ind);
    _dtde = ∂Θ∂E(sm, leaves, air, ind; δe = δe);

    #=
    if any(isnan, (_dade, _dtde))
        @info "Debugging" _dade _dtde;
        error("NaN detected in ∂g∂t calculation!");
    end;
    =#

    return max(-0.001, min(0.001, sm.K * (_dade - _dtde)))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT} = (
    return ∂g∂t(sm, leaves, air, sm.β.PARAM_Y, ind; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return max(-0.001, min(0.001, (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.τ))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return max(-0.001, min(0.001, (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.τ))
);


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

    ∂gₙ∂t(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT}

Return the marginal increase of stomatal conductance, given
- `lf` `Leaf`, `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions

"""
function ∂gₙ∂t end

∂gₙ∂t(lf::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = ∂gₙ∂t(lf.SM, lf, air);

∂gₙ∂t(sm::WangSM{FT}, lf::Leaves2D{FT}, air::AirLayer{FT}) where {FT} = max(-0.001, min(0.001, sm.K * (∂R∂E(lf, air) - ∂Θₙ∂E(lf, air))));


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""
This function updates stomatal conductance for H₂O and CO₂. Supported functionalities are
- Update conductance for H₂O prognostically
- Update conductance for CO₂ based on that for H₂O

"""
function stomatal_conductance! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-12: add new method to compute marginal stomatal conductance increase
#     2022-Jul-12: add new method to update marginal stomatal conductance for SPAC
#     2023-Sep-14: add root disconnection control
# To do
#     TODO: be careful with the β here (need to used the value stored in empirical SM)
#     TODO: use ∂gₙ∂t for nighttime conditions
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(spac::MultiLayerSPAC{FT}; β::FT = FT(1)) where {FT}

Update marginal stomatal conductance, given
- `spac` `MultiLayerSPAC` type struct
- `β` Tuning factor

"""
stomatal_conductance!(spac::MultiLayerSPAC{FT}; β::FT = FT(1)) where {FT} = (
    (; AIR, CANOPY, LEAVES, LEAVES_INDEX) = spac;

    # if lai = 0 or roots are not connected, do nothing
    # TODO: redo this later to foce dgdt to 0
    if CANOPY.lai == 0 || !spac._root_connection
        return nothing
    end;

    for _i in eachindex(LEAVES_INDEX)
        stomatal_conductance!(LEAVES[_i], AIR[LEAVES_INDEX[_i]]; β = β);
    end;

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    leaves.∂g∂t_shaded = ∂g∂t(leaves, air; β = β);
    for _i in eachindex(leaves.∂g∂t_sunlit)
        leaves.∂g∂t_sunlit[_i] = ∂g∂t(leaves, air, _i; β = β);
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add new method to update stomatal conductance for CO₂ based on that of H₂O
#     2022-Jul-07: add new method to update stomatal conductance prognostically
#     2022-Jul-12: move ∂g∂t to another method
#     2022-Jul-12: add method to update g for SPAC
#     2022-Jul-26: limit g in range after updating stomatal conductance
#     2023-Mar-11: do nothing if LAI == 0
#     2023-Mar-13: move some methods as stomatal_conductance_profile!
#     2023-Sep-14: add root disconnection control
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update stomatal conductance for H₂O based on computed ∂g∂t, given
- `spac` `MultiLayerSPAC` type struct
- `δt` Time step length `[s]`

"""
stomatal_conductance!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # if lai = 0 or roots are not connected, do nothing
    if CANOPY.lai == 0 || !spac._root_connection
        return nothing
    end;

    for _leaves in LEAVES
        stomatal_conductance!(_leaves, δt);
    end;

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}, δt::FT) where {FT} = (
    leaves.g_H₂O_s_shaded += leaves.∂g∂t_shaded * δt;
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        leaves.g_H₂O_s_sunlit[_i] += leaves.∂g∂t_sunlit[_i] * δt;
    end;
    limit_stomatal_conductance!(leaves);
    stomatal_conductance_profile!(leaves);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2023-Mar-13: add function to update stomatal conductance profile based on gs and gb
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance_profile!(spac::MultiLayerSPAC{FT}) where {FT}

Update stomatal conductance for CO₂ based on that for H₂O, given
- `spac` `MultiLayerSPAC` type struct

"""
function stomatal_conductance_profile! end

stomatal_conductance_profile!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    if CANOPY.lai == 0
        return nothing
    end;

    for _leaves in LEAVES
        stomatal_conductance_profile!(_leaves);
    end;

    return nothing
);

stomatal_conductance_profile!(leaves::Leaves2D{FT}) where {FT} = (
    leaves._g_CO₂_shaded = 1 / (1 / leaves.g_CO₂_b + FT(1.6) / leaves.g_H₂O_s_shaded);
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        leaves._g_CO₂_sunlit[_i] = 1 / (1 / leaves.g_CO₂_b + FT(1.6) / leaves.g_H₂O_s_sunlit[_i]);
    end;

    return nothing
);
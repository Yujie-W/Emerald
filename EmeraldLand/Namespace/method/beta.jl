# This file contains the beta functions

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for which parameter to tune
#     2022-Jun-30: add struct to tune G1
#     2022-Jun-30: add struct to base on Kleaf
#     2022-Jun-30: add struct to base on Ksoil
#     2022-Jun-30: add struct to base on Pleaf
#     2022-Jun-30: add struct to base on Psoil
#     2022-Jun-30: add struct to tune Vcmax
#     2022-Jun-30: add struct to base on Θ (SWC)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractBetaParameter`:
- `BetaParameterG1` PARAM_Y
- `BetaParameterKleaf` PARAM_X
- `BetaParameterKsoil` PARAM_X
- `BetaParameterPleaf` PARAM_X
- `BetaParameterPsoil` PARAM_X
- `BetaParameterVcmax` PARAM_Y
- `BetaParameterΘ` PARAM_X

"""
abstract type AbstractBetaParameter end;


""" Empty struct indicating that the function tunes G1 parameter of an empirical model """
struct BetaParameterG1 <: AbstractBetaParameter end;


""" Empty struct indicating that the beta function is based on Kleaf """
struct BetaParameterKleaf <: AbstractBetaParameter end;


""" Empty struct indicating that the beta function is based on Ksoil """
struct BetaParameterKsoil <: AbstractBetaParameter end;


""" Empty struct indicating that the beta function is based on Pleaf """
struct BetaParameterPleaf <: AbstractBetaParameter end;


""" Empty struct indicating that the beta function is based on Psoil """
struct BetaParameterPsoil <: AbstractBetaParameter end;


""" Empty struct indicating that the function tunes Vcmax for an empirical model """
struct BetaParameterVcmax <: AbstractBetaParameter end;


""" Empty struct indicating that the beta function is based on soil water content """
struct BetaParameterΘ <: AbstractBetaParameter end;

sync_state!(state_from::AbstractBetaParameter, state_to::AbstractBetaParameter) = nothing;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for modular beta function
#     2022-Jun-30: add more types to PARAM_X
#     2022-Jul-07: use BetaParameterKleaf as the default param_x
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to tune G1 or Vcmax based on leaf hydraulic conductance

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BetaFunction{FT<:AbstractFloat}
    # General model information
    "Function to turn variables to β tuning factor"
    FUNC::Function = (x -> x)
    "Input parameter to base on"
    PARAM_X::Union{BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ} = BetaParameterKleaf()
    "Target parameter to tune"
    PARAM_Y::Union{BetaParameterG1, BetaParameterVcmax} = BetaParameterG1()
end;

sync_state!(state_from::BetaFunction{FT}, state_to::BetaFunction{FT}) where {FT} = (
    state_to.FUNC = deepcopy(state_from.FUNC);
    sync_state!(state_from.PARAM_X, state_to.PARAM_X);
    sync_state!(state_from.PARAM_Y, state_to.PARAM_Y);

    return nothing
);

"""
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


"""
Struct to tune G1 or Vcmax based on leaf hydraulic conductance
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

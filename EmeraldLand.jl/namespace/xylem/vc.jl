# This file contains the definitions of the vulnerability curve models
# Because VC curve is a trait (a state), so there is no need to split the struct into state and auxilary parts.

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Feb-01: add abstract type for vulnerability curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractXylemVC:
- [`LogisticVC`](@ref)
- [`PowerVC`](@ref)
- [`WeibullVC`](@ref)
- [`ComplexVC`](@ref)

"""
abstract type AbstractXylemVC{FT<:AbstractFloat} end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add LogisticVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Modified logistic function for vulnerability curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LogisticVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Multiplier to exponential component"
    a::FT = 1
    "Multiplier to pressure `[MPa⁻¹]`"
    b::FT = 1
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add PowerVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Power function for vulnerability curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PowerVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Multiplier to power component `[MPa⁻ᵇ]`"
    a::FT = 1
    "Power to pressure"
    b::FT = 1
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add WeibullVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Weibull cumulative distribution function for vulnerability curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct WeibullVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Numerator in the exponential component `[MPa]`"
    b::FT = 2
    "Power to pressure component"
    c::FT = 5
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add ComplexVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

A complex struct for segmented vulnerability curve such as dual Weibull function

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ComplexVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Percentages of each VC component"
    fs::Vector{FT} = FT[0.5, 0.5]
    "Vector of vulnerability curve components"
    vcs::Union{Vector{LogisticVC{FT}}, Vector{PowerVC{FT}}, Vector{WeibullVC{FT}}, Vector{AbstractXylemVC{FT}}} = WeibullVC{FT}[WeibullVC{FT}(), WeibullVC{FT}(B = 3)]
end;

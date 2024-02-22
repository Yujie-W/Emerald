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
    A::FT = 1
    "Multiplier to pressure `[MPa⁻¹]`"
    B::FT = 1
end;

sync_state!(state_from::LogisticVC{FT}, state_to::LogisticVC{FT}) where {FT} = (
    state_to.A = state_from.A;
    state_to.B = state_from.B;

    return nothing
);


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
    A::FT = 1
    "Power to pressure"
    B::FT = 1
end;

sync_state!(state_from::PowerVC{FT}, state_to::PowerVC{FT}) where {FT} = (
    state_to.A = state_from.A;
    state_to.B = state_from.B;

    return nothing
);


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
    B::FT = 2
    "Power to pressure component"
    C::FT = 5
end;

sync_state!(state_from::WeibullVC{FT}, state_to::WeibullVC{FT}) where {FT} = (
    state_to.B = state_from.B;
    state_to.C = state_from.C;

    return nothing
);


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

sync_state!(state_from::ComplexVC{FT}, state_to::ComplexVC{FT}) where {FT} = (
    state_to.fs .= state_from.fs;
    for i in eachindex(state_from.vcs)
        sync_state!(state_from.vcs[i], state_to.vcs[i])
    end;

    return nothing
);

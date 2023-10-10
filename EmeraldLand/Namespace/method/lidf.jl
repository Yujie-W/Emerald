# This file contains structs for the leaf inclination angle distribution function (LIDF) algorithms

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-02: add abstract type for LIDF algorithms
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractLIDFAlgorithm:
- [`BetaLIDF`](@ref)
- [`VerhoefLIDF`](@ref)

"""
abstract type AbstractLIDFAlgorithm{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2023-May-22: add beta function method
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Beta LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BetaLIDF{FT<:AbstractFloat} <: AbstractLIDFAlgorithm{FT}
    # General model information
    "Leaf inclination angle distribution function parameter a"
    A::FT = 1
    "Leaf inclination angle distribution function parameter b"
    B::FT = 1
end;


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to HyperspectralMLCanopy
#     2022-Jun-02: abstractize LIDF as a field
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VerhoefLIDF{FT<:AbstractFloat} <: AbstractLIDFAlgorithm{FT}
    # General model information
    "Leaf inclination angle distribution function parameter a"
    A::FT = 0
    "Leaf inclination angle distribution function parameter b"
    B::FT = 0
end;

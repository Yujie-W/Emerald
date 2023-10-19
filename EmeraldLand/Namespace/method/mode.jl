# This file contains the photosynthesis mode

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-14: add abstract mode type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractPhotosynthesisMode:
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)

"""
abstract type AbstractPhotosynthesisMode end;


""" An empty structure to signal the function to calculate photosynthetic rates based on leaf diffusive conductance to CO₂ """
struct GCO₂Mode <: AbstractPhotosynthesisMode end;


""" An empty structure to signal the function to calculate photosynthetic rates based on CO₂ partial pressure """
struct PCO₂Mode <: AbstractPhotosynthesisMode end;

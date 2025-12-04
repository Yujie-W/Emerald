#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2024-02-05: add abstract type for generalized flux tower
#     2024-02-05: add US_NR1 struct as a type of abstract flux tower
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Abstract flux tower, each flux tower may have different data structure

"""
abstract type AbstractFluxTower end;

struct US_NR1 <: AbstractFluxTower end;

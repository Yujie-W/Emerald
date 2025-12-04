#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2024-Sep-07: add a structure to save canopy clumping index related parameters
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save canopy clumping index

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ClumpingIndex{FT}
    "Clumping index when zenith angle is 0 degree"
    ci_0::FT = 1.0
    "Tuning factor to account for the impact of zenith angle on clumping index"
    ci_1::FT = 0.0
end

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Aug-29: add carbon pool
#     2024-Aug-30: add max and min threshold for carbon pool (portion exceed max must be used for new growth; when below min, no new growth or recovery)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for carbon pool (one)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CarbonPoolWholePlant{FT}
    "Carbon pool `[mol]`"
    c_pool::FT = 3 * 80 * 0.02 * 10000 / 30
    "Maximum threshold of pool, extra must be used for new growth (need to scale with biomass?) `[mol]`"
    c_pool_max::FT = 3 * 80 * 0.02 * 10000 / 30 * 10
    "Minimal threshold for new growth (need to scale with leaf biomass) `[mol]`"
    c_pool_min::FT = 3 * 80 * 0.02 * 10000 / 30
end;

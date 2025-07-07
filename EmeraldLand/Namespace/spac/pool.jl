#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Aug-29: add carbon pool
#     2024-Aug-30: add max and min threshold for carbon pool (portion exceed max must be used for new growth; when below min, no new growth or recovery)
#     2024-Sep-03: remove field of min threshold for carbon pool
#     2024-Sep-04: add field of min threshold back so that carbon pool does not go to zero immediately when leaf regrow flag is false
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for carbon pool (one)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CarbonPoolWholePlant{FT}
    "Carbon pool (default is one time the LAI investment at max LAI) `[mol]`"
    c_pool::FT = 3 * 80 * 0.02 * 10000 / 30
    "Maximum threshold of pool, extra must be used for new growth (default is 2 times the LAI + min; need to scale with biomass?) `[mol]`"
    c_pool_max::FT = 2.5 * c_pool
    "Minimum threshold of pool, new LAI growth allower unless _leaf_regrow flag is true `[mol]`"
    c_pool_min::FT = 0.5 * c_pool
end;

sync_state!(state_from::CarbonPoolWholePlant{FT}, state_to::CarbonPoolWholePlant{FT}) where {FT} = (
    state_to.c_pool = state_from.c_pool;
    state_to.c_pool_max = state_from.c_pool_max;
    state_to.c_pool_min = state_from.c_pool_min;

    return nothing
);

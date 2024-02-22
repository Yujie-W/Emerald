# This file contains the state and auxil structs for the SPAC memory (auxil vars as cache variables)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct PlantMemoryState
#     2024-Feb-22: remove state and auxil from memory struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state and state variables of the SPAC memory

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PlantMemory{FT}
    "Chlorophyll content used for last time step"
    chl::FT = -1
    "LAI used for last time step"
    lai::FT = -1
    "Memory about the historical temperature"
    t_history::Vector{FT} = ones(FT,10) .* FT(NaN)
    "Vcmax25 used for last time step"
    vcmax25::FT = -1
end;

sync_state!(state_from::PlantMemory{FT}, state_to::PlantMemory{FT}) where {FT} = (
    state_to.chl = state_from.chl;
    state_to.lai = state_from.lai;
    state_to.t_history .= state_from.t_history;
    state_to.vcmax25 = state_from.vcmax25;

    return nothing
);

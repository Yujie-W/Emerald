# This file contains the state and auxil structs for the SPAC memory (auxil vars as cache variables)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct PlantMemoryState
#     2024-Feb-22: remove state and auxil from memory struct
#     2024-Jul-24: remove lai, ci, vcmax, and cab from memory (use traits instead)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state and state variables of the SPAC memory

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PlantMemory{FT}
    "Memory about the historical temperature"
    t_history::Vector{FT}
    "Index of the current time step"
    i_history::Int = 1
end;

PlantMemory(config::SPACConfiguration{FT}) where {FT} = PlantMemory{FT}(t_history = ones(FT, config.T_CLM_N) .* T₂₅());

sync_state!(state_from::PlantMemory{FT}, state_to::PlantMemory{FT}) where {FT} = (
    state_to.t_history .= state_from.t_history;
    state_to.i_history = state_from.i_history;

    return nothing
);

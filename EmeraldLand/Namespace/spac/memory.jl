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
    t_history::Vector{FT} = FT[298]
end;

sync_state!(state_from::PlantMemory{FT}, state_to::PlantMemory{FT}) where {FT} = (
    state_to.t_history = deepcopy(state_from.t_history);

    return nothing
);

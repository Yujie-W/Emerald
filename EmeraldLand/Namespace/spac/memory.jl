# This file contains the state and auxil structs for the SPAC memory (auxil vars as cache variables)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct PlantMemoryState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state and state variables of the SPAC memory

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PlantMemoryState{FT}
    "Memory about the historical temperature"
    t_history::Vector{FT} = ones(FT,10) .* FT(NaN)
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct PlantMemoryAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the auxiliary and auxiliary variables of the SPAC memory

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PlantMemoryAuxil{FT}
    "Chlorophyll content used for last time step"
    chl::FT = -1
    "LAI used for last time step"
    lai::FT = -1
    "Vcmax25 used for last time step"
    vcmax25::FT = -1
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct BulkSPACMemory
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state and auxiliary variables of the SPAC memory

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PlantMemory{FT}
    "State variables"
    state::PlantMemoryState{FT} = PlantMemoryState{FT}()
    "Auxiliary variables"
    auxil::PlantMemoryAuxil{FT} = PlantMemoryAuxil{FT}()
end;

# This file contains the structs for the plant in SPAC

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct Plant
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for plant in SPAC

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Plant{FT}
    "Plant root depth and height information `[m]`"
    zs::Vector{FT}

    "Root hydraulic system"
    roots::Vector{Root{FT}}
    "Corresponding soil layer per root layer"
    roots_index::Vector{Int}
    "Trunk hydraulic system"
    trunk::Stem{FT}
    "Root-trunk junction capacitor used for roots flow calculations"
    junction::JunctionCapacitor{FT} = JunctionCapacitor{FT}()
    "Branch hydraulic system"
    branches::Vector{Stem{FT}}
    "Leaf per layer"
    leaves::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    leaves_index::Vector{Int}
    "Memory cache"
    memory::PlantMemory{FT} = PlantMemory{FT}()

    # Cache variables
    "Whether there is any root connected to soil"
    _root_connection::Bool = true
end;

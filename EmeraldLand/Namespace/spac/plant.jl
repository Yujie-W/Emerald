# This file contains the structs for the plant in SPAC

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct Plant
#     2024-Jul-24: add leaf shedded flag cache
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
    "Root-trunk junction capacitor used for roots flow calculations"
    junction::JunctionCapacitor{FT} = JunctionCapacitor{FT}()
    "Trunk hydraulic system"
    trunk::Stem{FT}
    "Branch hydraulic system"
    branches::Vector{Stem{FT}}
    "Leaf per layer"
    leaves::Vector{CanopyLayer{FT}}
    "Corresponding air layer per canopy layer"
    leaves_index::Vector{Int}
    "Memory cache"
    memory::PlantMemory{FT}

    # Cache variables
    "Whether leaves are shedded"
    _leaf_shedded::Bool = false
    "Whether there is any root connected to soil"
    _root_connection::Bool = true
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-22: add struct PlantStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save plant states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct PlantStates{FT<:AbstractFloat}
    "Plant root depth and height information `[m]`"
    zs::Vector{FT}
    "Root hydraulic system"
    roots::Vector{RootStates{FT}}
    "Corresponding soil layer per root layer"
    roots_index::Vector{Int}
    "Root-trunk junction capacitor used for roots flow calculations"
    junction::JunctionCapacitorState{FT}
    "Trunk hydraulic system"
    trunk::StemStates{FT}
    "Branch hydraulic system"
    branches::Vector{StemStates{FT}}
    "CanopyLayer or Leaf per layer"
    leaves::Vector{CanopyLayerStates{FT}}
    "Corresponding air layer per canopy layer"
    leaves_index::Vector{Int}
    "Memory cache"
    memory::PlantMemory{FT}
end;

PlantStates(plant::Plant{FT}) where {FT} = PlantStates{FT}(
            plant.zs,
            [RootStates(root) for root in plant.roots],
            plant.roots_index,
            plant.junction.state,
            StemStates(plant.trunk),
            [StemStates(branch) for branch in plant.branches],
            [CanopyLayerStates(l) for l in plant.leaves],
            plant.leaves_index,
            plant.memory
);

sync_state!(plant::Plant{FT}, states::PlantStates{FT}) where {FT} = (
    states.zs .= plant.zs;
    for i in 1:length(plant.roots)
        sync_state!(plant.roots[i], states.roots[i]);
    end;
    states.roots_index .= plant.roots_index;
    sync_state!(plant.junction.state, states.junction);
    sync_state!(plant.trunk, states.trunk);
    for i in 1:length(plant.branches)
        sync_state!(plant.branches[i], states.branches[i]);
    end;
    for i in 1:length(plant.leaves)
        sync_state!(plant.leaves[i], states.leaves[i]);
    end;
    states.leaves_index .= plant.leaves_index;
    sync_state!(plant.memory, states.memory);

    return nothing
);

sync_state!(states::PlantStates{FT}, plant::Plant{FT}) where {FT} = (
    plant.zs .= states.zs;
    for i in 1:length(plant.roots)
        sync_state!(states.roots[i], plant.roots[i]);
    end;
    plant.roots_index .= states.roots_index;
    sync_state!(states.junction, plant.junction.state);
    sync_state!(states.trunk, plant.trunk);
    for i in 1:length(plant.branches)
        sync_state!(states.branches[i], plant.branches[i]);
    end;
    for i in 1:length(plant.leaves)
        sync_state!(states.leaves[i], plant.leaves[i]);
    end;
    plant.leaves_index .= states.leaves_index;
    sync_state!(states.memory, plant.memory);

    return nothing
);

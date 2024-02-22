# This file contains the fields within a leaf

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: add Leaf struct
#     2023-Oct-03: add fields photosystem, flux
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the fields within a leaf

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaf{FT}
    "Leaf biophysics struct"
    bio::LeafBio{FT}
    "Extraxylary capacitor struct"
    capacitor::ExtraXylemCapacitor{FT} = ExtraXylemCapacitor{FT}()
    "Leaf energy struct"
    energy::LeafEnergy{FT} = LeafEnergy{FT}()
    "Leaf flux struct"
    flux::LeafFlux{FT}
    "Photosynthesis system struct"
    photosystem::LeafPhotosystem{FT} = LeafPhotosystem{FT}()
    "Leaf xylem hydraulics struct"
    xylem::XylemHydraulics{FT}
end;


"""

    Leaf(config::SPACConfiguration{FT}) where {FT}

Return the leaf struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Leaf(config::SPACConfiguration{FT}) where {FT} = Leaf{FT}(bio = LeafBio(config), flux = LeafFlux(config), xylem = XylemHydraulics(config));


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-22: add struct LeafStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save leaf states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct LeafStates{FT<:AbstractFloat}
    "Leaf biophysics state"
    bio::LeafBioState{FT}
    "Extraxylary capacitor state"
    capacitor::ExtraXylemCapacitorState{FT}
    "Leaf energy state"
    energy::LeafEnergyState{FT}
    "Leaf flux state"
    flux::LeafFluxState{FT}
    "Photosynthesis system state"
    photosystem::Union{C3CytoState{FT}, C3VJPState{FT}, C4VJPState{FT}}
    "Leaf xylem hydraulics state"
    xylem::XylemHydraulicsState{FT}
end;

LeafStates(leaf::Leaf{FT}) where {FT} = LeafStates{FT}(leaf.bio.state, leaf.capacitor.state, leaf.energy.state, leaf.flux.state, leaf.photosystem.state, leaf.xylem.state);

sync_state!(leaf::Leaf{FT}, states::LeafStates{FT}) where {FT} = (
    sync_state!(leaf.bio.state, states.bio);
    sync_state!(leaf.capacitor.state, states.capacitor);
    sync_state!(leaf.energy.state, states.energy);
    sync_state!(leaf.flux.state, states.flux);
    sync_state!(leaf.photosystem.state, states.photosystem);
    sync_state!(leaf.xylem.state, states.xylem);

    return nothing
);

sync_state!(states::LeafStates{FT}, leaf::Leaf{FT}) where {FT} = (
    sync_state!(states.bio, leaf.bio.state);
    sync_state!(states.capacitor, leaf.capacitor.state);
    sync_state!(states.energy, leaf.energy.state);
    sync_state!(states.flux, leaf.flux.state);
    sync_state!(states.photosystem, leaf.photosystem.state);
    sync_state!(states.xylem, leaf.xylem.state);

    return nothing
);

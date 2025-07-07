# This file contains the fields within a canopy layer

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-25: add struct CanopyLayer
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the fields within a canopy layer

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyLayer{FT}
    "CanopyLayer biophysics struct"
    bio::LeafBio{FT}
    "Extraxylary capacitor struct"
    capacitor::ExtraXylemCapacitor{FT} = ExtraXylemCapacitor{FT}()
    "CanopyLayer energy struct"
    energy::LeafEnergy{FT} = LeafEnergy{FT}()
    "CanopyLayer flux struct"
    flux::CanopyLayerFlux{FT}
    "Photosynthesis system struct"
    photosystem::CanopyLayerPhotosystem{FT}
    "CanopyLayer xylem hydraulics struct"
    xylem::XylemHydraulics{FT}
end;


"""

    CanopyLayer(config::SPACConfiguration{FT}) where {FT}

Return the CanopyLayer struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
CanopyLayer(config::SPACConfiguration{FT}) where {FT} = (
    clayer = CanopyLayer{FT}(bio = LeafBio(config), flux = CanopyLayerFlux(config), photosystem = CanopyLayerPhotosystem(config), xylem = XylemHydraulics(config));
    clayer.xylem.trait.cp = 1780;
    clayer.xylem.trait.k_max = 0.04;

    return clayer
);

kill_plant!(st::CanopyLayer{FT}) where {FT} = (
    kill_plant!(st.bio);
    kill_plant!(st.capacitor);
    kill_plant!(st.energy);
    kill_plant!(st.flux);
    kill_plant!(st.photosystem);
    kill_plant!(st.xylem);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-25: add struct CanopyLayerStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save canopy layer states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct CanopyLayerStates{FT<:AbstractFloat}
    "CanopyLayer biophysics state"
    bio::LeafBioState{FT}
    "Extraxylary capacitor state"
    capacitor::ExtraXylemCapacitorState{FT}
    "CanopyLayer energy state"
    energy::LeafEnergyState{FT}
    "CanopyLayer flux state"
    flux::CanopyLayerFluxState{FT}
    "Photosynthesis system state"
    photosystem::Union{C3State{FT}, C4State{FT}}
    "CanopyLayer xylem hydraulics state"
    xylem::XylemHydraulicsState{FT}
end;

CanopyLayerStates(clayer::CanopyLayer{FT}) where {FT} = CanopyLayerStates{FT}(
            deepcopy(clayer.bio.state),
            deepcopy(clayer.capacitor.state),
            deepcopy(clayer.energy.state),
            deepcopy(clayer.flux.state),
            deepcopy(clayer.photosystem.state),
            deepcopy(clayer.xylem.state));

sync_state!(clayer::CanopyLayer{FT}, states::CanopyLayerStates{FT}) where {FT} = (
    sync_state!(clayer.bio.state, states.bio);
    sync_state!(clayer.capacitor.state, states.capacitor);
    sync_state!(clayer.energy.state, states.energy);
    sync_state!(clayer.flux.state, states.flux);
    sync_state!(clayer.photosystem.state, states.photosystem);
    sync_state!(clayer.xylem.state, states.xylem);

    return nothing
);

sync_state!(states::CanopyLayerStates{FT}, clayer::CanopyLayer{FT}) where {FT} = (
    sync_state!(states.bio, clayer.bio.state);
    sync_state!(states.capacitor, clayer.capacitor.state);
    sync_state!(states.energy, clayer.energy.state);
    sync_state!(states.flux, clayer.flux.state);
    sync_state!(states.photosystem, clayer.photosystem.state);
    sync_state!(states.xylem, clayer.xylem.state);

    return nothing
);

# This file contains stem state and auxil structs

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-25: add new Stem struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save stem energy and hydraulics variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Stem{FT}
    "Stem energy struct"
    energy::XylemEnergy{FT} = XylemEnergy{FT}()
    "Stem hydraulics struct"
    xylem::XylemHydraulics{FT}
end;

Stem(config::SPACConfiguration{FT}) where {FT} = Stem{FT}(xylem = XylemHydraulics(config));


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-22: add struct StemStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save stem states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct StemStates{FT<:AbstractFloat}
    "Stem energy state"
    energy::XylemEnergyState{FT}
    "Stem hydraulics state"
    xylem::XylemHydraulicsState{FT}
end;

StemStates(stem::Stem{FT}) where {FT} = StemStates{FT}(stem.energy.state, stem.xylem.state);

sync_state!(stem::Stem{FT}, states::StemStates{FT}) where {FT} = (
    sync_state!(stem.energy.state, states.energy);
    sync_state!(stem.xylem.state, states.xylem);

    return nothing
);

sync_state!(states::StemStates{FT}, stem::Stem{FT}) where {FT} = (
    sync_state!(states.energy, stem.energy.state);
    sync_state!(states.xylem, stem.xylem.state);

    return nothing
);

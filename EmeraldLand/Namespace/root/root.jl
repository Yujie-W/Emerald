# This file contains the root struct as a combination of the xylem and energy structs
# Note here that the energy struct needs to be initialized with the heat capacity of root and water, water content, and temperature...

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add Root struct with energy and xylem fields
#     2023-Sep-23: add constructor for Root struct and initialize the energy state of the root
#     2023-Sep-23: add field rhizosphere
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Root{FT}
    "Root energy struct"
    energy::XylemEnergy{FT} = XylemEnergy{FT}()
    "Rhizosphere struct"
    rhizosphere::Rhizosphere{FT} = Rhizosphere{FT}()
    "Root xylem struct"
    xylem::XylemHydraulics{FT}
end;

Root(config::SPACConfiguration{FT}) where {FT} = Root{FT}(xylem = XylemHydraulics(config));


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-22: add struct RootStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save root states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct RootStates{FT<:AbstractFloat}
    "Root energy state"
    energy::XylemEnergyState{FT}
    "Rhizosphere state"
    rhizosphere::RhizosphereState{FT}
    "Root xylem state"
    xylem::XylemHydraulicsState{FT}
end;

RootStates(root::Root{FT}) where {FT} = RootStates{FT}(root.energy.state, root.rhizosphere.state, root.xylem.state);

sync_state!(root::Root{FT}, states::RootStates{FT}) where {FT} = (
    sync_state!(root.energy.state, states.energy);
    sync_state!(root.rhizosphere.state, states.rhizosphere);
    sync_state!(root.xylem.state, states.xylem);

    return nothing
);

sync_state!(states::RootStates{FT}, root::Root{FT}) where {FT} = (
    sync_state!(states.energy, root.energy.state);
    sync_state!(states.rhizosphere, root.rhizosphere.state);
    sync_state!(states.xylem, root.xylem.state);

    return nothing
);

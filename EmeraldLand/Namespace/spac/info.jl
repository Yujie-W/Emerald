# This file contains structs to store general SPAC information such as geometry and site information (e.g. latitude, longitude, elevation)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct SPACInfo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for general SPAC information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACInfo{FT}
    # Geometry information
    "Air boundaries `[m]`"
    z_air::Vector{FT}
    "Soil boundaries `[m]`"
    z_soil::Vector{FT}

    # Geographical information
    "Elevation"
    elev::FT = 100
    "Latitude"
    lat::FT
    "Longitude"
    lon::FT
end;

sync_state!(state_from::SPACInfo{FT}, state_to::SPACInfo{FT}) where {FT} = (
    state_to.z_air .= state_from.z_air;
    state_to.z_soil .= state_from.z_soil;
    state_to.elev = state_from.elev;
    state_to.lat = state_from.lat;
    state_to.lon = state_from.lon;

    return nothing
);

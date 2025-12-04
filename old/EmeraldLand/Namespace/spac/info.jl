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

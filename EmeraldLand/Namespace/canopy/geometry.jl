#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: add extra fields haa, hsa, saa, and vaa, and remove field raa (relative azimuth angle, will be computed based on saa and vaa)
#     2022-Jul-20: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores sun sensor geometry information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunSensorGeometryToDelete{FT<:AbstractFloat}
    # General site information
    "Hill facing azimuth angle `[°]`, 0 for North, 180 for south"
    haa::FT = 0
    "Hill slope angle `[°]`"
    hsa::FT = 0
end

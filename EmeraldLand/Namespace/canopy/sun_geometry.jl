# This file contains structs for sun geometry (apart from sensor geometry)

#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SunGeometryState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state variables of the sun geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunGeometryState{FT}
    "Solar azimuth angle `[°]`, a function of time"
    saa::FT = 180
    "Solar zenith angle `[°]`, a function of lat and time"
    sza::FT = 30
end;


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SunGeometryAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the auxiliary variables of the sun geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunGeometryAuxil{FT}
end;


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SunGeometry
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the sun geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunGeometry{FT}
    "State variables"
    state::SunGeometryState{FT} = SunGeometryState{FT}()
    "Auxiliary variables"
    auxil::SunGeometryAuxil{FT} = SunGeometryAuxil{FT}()
end;

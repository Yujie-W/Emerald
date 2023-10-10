# This file contains structs for sensor geometry (apart from sun geometry)

#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SensorGeometryState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state variables of the sensor geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SensorGeometryState{FT}
    "Viewing azimuth angle `[°]`, a function of time"
    vaa::FT = 180
    "Viewing zenith angle `[°]`, a function of lat and time"
    vza::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SensorGeometryAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the auxiliary variables of the sensor geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SensorGeometryAuxil{FT}
end;


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SensorGeometry
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the sensor geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SensorGeometry{FT}
    "State variables"
    state::SensorGeometryState{FT} = SensorGeometryState{FT}()
    "Auxiliary variables"
    auxil::SensorGeometryAuxil{FT} = SensorGeometryAuxil{FT}()
end;

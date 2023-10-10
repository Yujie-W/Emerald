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
    # Scattering coefficients
    "Backward diffuse->observer scatter weight"
    dob::FT = 0
    "Forward diffuse->observer scatter weight"
    dof::FT = 0
    "Backward direct->observer scatter weight"
    sob::FT = 0
    "Forward direct->observer scatter weight"
    sof::FT = 0

    # Extinction coefficient related
    "Observer direction beam extinction coefficient weight (diffuse)"
    ko::FT = 0
    "Probability of directly viewing a leaf in observer direction at different layer boundaries"
    po::Vector{FT}
    "Bi-directional probability of directly viewing a leaf at different layer boundaries (solar->canopy->observer)"
    pso::Vector{FT}
end;

SensorGeometryAuxil(config::SPACConfiguration{FT}) where {FT} = SensorGeometryAuxil{FT}(
            po  = Vector{FT}(undef, config.DIM_LAYER + 1),
            pso = Vector{FT}(undef, config.DIM_LAYER + 1)
);


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
    auxil::SensorGeometryAuxil{FT}
end;

SensorGeometry(config::SPACConfiguration{FT}) where {FT} = SensorGeometry{FT}(auxil = SensorGeometryAuxil(config));

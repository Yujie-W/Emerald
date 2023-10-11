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

    # Extinction coefficient related (for different inclination angles)
    "cos(inclination) * cos(vza) at different inclination angles"
    Co_incl::Vector{FT}
    "sin(inclination) * sin(vza) at different inclination angles"
    So_incl::Vector{FT}
    "Outgoing beam extinction coefficient weights at different inclination angles"
    ko_incl::Vector{FT}
    "Backward scattering coefficients at different inclination angles"
    sb_incl::Vector{FT}
    "Forward scattering coefficients at different inclination angles"
    sf_incl::Vector{FT}
    "Co >= So ? FT(π) : acos(-Co/So)"
    βo_incl::Vector{FT}

    # Matrix used for radiation to sensor
    "Conversion factor fo for angle towards observer at different inclination and azimuth angles"
    fo::Matrix{FT}
    "Absolute value of fo"
    fo_abs::Matrix{FT}
    "fo * cos² Θ_INCL"
    fo_cos²_incl::Matrix{FT}
    "fo * fs"
    fo_fs::Matrix{FT}
    "Absolute value of fo * fs"
    fo_fs_abs::Matrix{FT}

    # Scattering coefficients per leaf area
    "Backward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    dob_leaf::Matrix{FT}
    "Forward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    dof_leaf::Matrix{FT}
    "Bidirectional from solar to observer scattering coefficient at different layers and wavelength bins"
    so_leaf::Matrix{FT}
end;

SensorGeometryAuxil(config::SPACConfiguration{FT}) where {FT} = SensorGeometryAuxil{FT}(
            po           = zeros(FT, config.DIM_LAYER + 1),
            pso          = zeros(FT, config.DIM_LAYER + 1),
            Co_incl      = zeros(FT, config.DIM_INCL),
            So_incl      = zeros(FT, config.DIM_INCL),
            ko_incl      = zeros(FT, config.DIM_INCL),
            sb_incl      = zeros(FT, config.DIM_INCL),
            sf_incl      = zeros(FT, config.DIM_INCL),
            βo_incl      = zeros(FT, config.DIM_INCL),
            fo           = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_abs       = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_cos²_incl = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_fs        = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_fs_abs    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            dob_leaf     = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            dof_leaf     = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            so_leaf      = zeros(FT, config.DIM_WL, config.DIM_LAYER),
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

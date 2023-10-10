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
    # Scattering coefficients
    "Backward diffuse->diffuse scatter weight"
    ddb::FT = 0
    "Forward diffuse->diffuse scatter weight"
    ddf::FT = 0
    "Backward direct->diffuse scatter weight"
    sdb::FT = 0
    "Forward direct->diffuse scatter weight"
    sdf::FT = 0

    # Extinction coefficient related
    "Solar direction beam extinction coefficient weight (direct)"
    ks::FT = 0
    "Probability of directly viewing a leaf in solar direction at different layers"
    p_sunlit::Vector{FT}
    "Probability of directly viewing a leaf in solar direction at different layer boundaries"
    ps::Vector{FT}

    # Extinction coefficient related (for different inclination angles)
    "cos(inclination) * cos(sza) at different inclination angles"
    Cs_incl::Vector{FT}
    "sin(inclination) * sin(sza) at different inclination angles"
    Ss_incl::Vector{FT}
    "Solar beam extinction coefficient weights at different inclination angles"
    ks_incl::Vector{FT}
    "Cs >= Ss ? FT(π) : acos(-Cs/Ss)"
    βs_incl::Vector{FT}

    # Matrix used for solar radiation
    "Conversion factor fs for angles from solar at different inclination and azimuth angles"
    fs::Matrix{FT}
    "Absolute value of fs"
    fs_abs::Matrix{FT}
    "fs * cos Θ_INCL"
    fs_cos²_incl::Matrix{FT}
end;

SunGeometryAuxil(config::SPACConfiguration{FT}) where {FT} = SunGeometryAuxil{FT}(
            p_sunlit     = zeros(FT, config.DIM_LAYER),
            ps           = zeros(FT, config.DIM_LAYER + 1),
            Cs_incl      = zeros(FT, config.DIM_INCL),
            Ss_incl      = zeros(FT, config.DIM_INCL),
            ks_incl      = zeros(FT, config.DIM_INCL),
            βs_incl      = zeros(FT, config.DIM_INCL),
            fs           = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fs_abs       = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fs_cos²_incl = zeros(FT, config.DIM_INCL, config.DIM_AZI),
);


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
    auxil::SunGeometryAuxil{FT}
end;

SunGeometry(config::SPACConfiguration{FT}) where {FT} = SunGeometry{FT}(auxil = SunGeometryAuxil(config));

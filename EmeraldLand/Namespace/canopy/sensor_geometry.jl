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
#     2024-Feb-25: add struct SensorGeometrySDAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state-dependetn auxiliary variables of the sensor geometry

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SensorGeometrySDAuxil{FT}
    # Scattering coefficients
    "Backward diffuse->observer scatter weight"
    dob::FT = 0
    "Forward diffuse->observer scatter weight"
    dof::FT = 0
    "Backward direct->observer scatter weight"
    sob::FT = 0
    "Forward direct->observer scatter weight"
    sof::FT = 0

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

    # Extinction coefficient related
    "Observer direction beam extinction coefficient weight (diffuse)"
    ko::FT = 0
    "Probability of directly viewing a leaf in observer direction at different layer boundaries"
    p_sensor::Vector{FT}
    "Probability of directly viewing soil in observer direction at different layer boundaries"
    p_sensor_soil::FT = 0
    "Bi-directional probability of directly viewing a leaf at different layer boundaries (solar->canopy->observer)"
    p_sun_sensor::Vector{FT}

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
end;

SensorGeometrySDAuxil(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = SensorGeometrySDAuxil{FT}(
            Co_incl      = zeros(FT, config.DIM_INCL),
            So_incl      = zeros(FT, config.DIM_INCL),
            ko_incl      = zeros(FT, config.DIM_INCL),
            sb_incl      = zeros(FT, config.DIM_INCL),
            sf_incl      = zeros(FT, config.DIM_INCL),
            βo_incl      = zeros(FT, config.DIM_INCL),
            p_sensor     = zeros(FT, n_layer),
            p_sun_sensor = zeros(FT, n_layer),
            fo           = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_abs       = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_cos²_incl = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_fs        = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fo_fs_abs    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
);


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SensorGeometryAuxil
#     2023-Oct-18: add fields dob_stem, dof_stem, and so_stem
#     2024-Jul-27: use bined PPAR to speed up (moved sif yield here from leaf)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the auxiliary variables of the sensor geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SensorGeometryAuxil{FT}
    # Scattering coefficients per leaf area or stem area
    "Backward scattering coefficient for diffuse->observer at different layers and wavelength bins of leaf"
    dob_leaf::Matrix{FT}
    "Forward scattering coefficient for diffuse->observer at different layers and wavelength bins of leaf"
    dof_leaf::Matrix{FT}
    "Bidirectional from solar to observer scattering coefficient at different layers and wavelength bins of leaf"
    so_leaf::Matrix{FT}
    "Backward scattering coefficient for diffuse->observer at different layers and wavelength bins of stem"
    dob_stem::Matrix{FT}
    "Forward scattering coefficient for diffuse->observer at different layers and wavelength bins of stem"
    dof_stem::Matrix{FT}
    "Bidirectional from solar to observer scattering coefficient at different layers and wavelength bins of stem"
    so_stem::Matrix{FT}

    # Canopy reflection
    "Radiation towards the viewing direction per layer (including soil) `[mW m⁻² nm⁻¹]`"
    e_sensor_layer::Matrix{FT}
    "Total radiation towards the viewing direction `[mW m⁻² nm⁻¹]`"
    e_sensor::Vector{FT}
    "Reflectance of the soil-canopy at the sensor direction"
    reflectance::Vector{FT}

    # Canopy fluorescence
    "SIF that can be observed from scattered SIF"
    sif_scattered::Matrix{FT}
    "SIF that can be observed from the shaded part of the canopy"
    sif_shaded::Matrix{FT}
    "SIF that can be observed from the sunlit part of the canopy"
    sif_sunlit::Matrix{FT}
    "SIF at sensor direction contributed by shaded leaves fo (not hotspot)"
    sif_obs_shaded::Vector{FT}
    "SIF at sensor direction contributed from soil"
    sif_obs_soil::Vector{FT}
    "SIF at sensor direction contributed by sunlit leaves fo_fs (hotspot)"
    sif_obs_sunlit::Vector{FT}
    "SIF at sensor direction contributed by scattered SIF"
    sif_obs_scattered::Vector{FT}
    "Total SIF at sensor direction"
    sif_obs::Vector{FT}
    "Shaded SIF yield of all ther layers"
    ϕ_f_shaded::Vector{FT}
    "Sunlit SIF yield of all ther layers"
    ϕ_f_sunlit::Vector{Matrix{FT}}
end;

SensorGeometryAuxil(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = SensorGeometryAuxil{FT}(
            dob_leaf          = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            dof_leaf          = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            so_leaf           = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            dob_stem          = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            dof_stem          = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            so_stem           = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            e_sensor_layer    = zeros(FT, length(config.SPECTRA.Λ), n_layer + 1),
            e_sensor          = zeros(FT, length(config.SPECTRA.Λ)),
            reflectance       = zeros(FT, length(config.SPECTRA.Λ)),
            sif_scattered     = zeros(FT, length(config.SPECTRA.IΛ_SIF), n_layer),
            sif_shaded        = zeros(FT, length(config.SPECTRA.IΛ_SIF), n_layer),
            sif_sunlit        = zeros(FT, length(config.SPECTRA.IΛ_SIF), n_layer),
            sif_obs_scattered = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
            sif_obs_shaded    = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
            sif_obs_sunlit    = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
            sif_obs_soil      = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
            sif_obs           = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
            ϕ_f_shaded        = zeros(FT, n_layer),
            ϕ_f_sunlit        = Matrix{FT}[ zeros(FT, config.DIM_INCL, config.DIM_AZI) for i in 1:n_layer],
);


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SensorGeometry
#     2023-Oct-18: add field s_aux
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the sensor geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SensorGeometry{FT}
    "State variables that may evolve with time"
    state::SensorGeometryState{FT} = SensorGeometryState{FT}()
    "State-dependent variables"
    s_aux::SensorGeometrySDAuxil{FT}
    "Auxiliary variables"
    auxil::SensorGeometryAuxil{FT}
end;

SensorGeometry(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = SensorGeometry{FT}(
            s_aux = SensorGeometrySDAuxil(config, n_layer),
            auxil = SensorGeometryAuxil(config, n_layer)
);

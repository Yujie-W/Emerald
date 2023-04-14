#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-13: add state struct to save SPAC configurations
#     2023-Apr-13: move Φ_PHOTON, RAD_SW_REF from MultiLayerSPAC
#     2023-Apr-13: move APAR_CAR from leaf structs
#     2023-Apr-13: move LHA and WLSET from HyperspectralMLCanopy
#     2023-Apr-13: add fields DIM_* and STEADY_STATE_HS
#     2023-Apr-13: move MAT_ρ from HyperspectralSoilAlbedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Global configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACConfiguration{FT}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # Constants
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT} = FT[700, 2500]
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT} = FT[400, 750]
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT} = FT[640, 850]
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT} = FT[400, 750]

    # Dimensions
    DIMS::SPACDimension = SPACDimension(DATASET, WL_NIR, WL_PAR, WL_SIF, WL_SIFE)

    # Wavelength information
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT} = WaveLengthSet{FT,DIMS}(DATASET, WL_NIR, WL_PAR, WL_SIF, WL_SIFE)

    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to use steady state plant hydraulic system"
    STEADY_STATE_HS::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Vectors
    "A matrix of characteristic curves"
    MAT_ρ::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]

    # Embedded structures
    "Hyperspectral absorption features of different leaf components"
    LHA::HyperspectralAbsorption{FT} = HyperspectralAbsorption{FT,DIMS}(DATASET)
    "Downwelling shortwave radiation reference spectrum"
    RAD_SW_REF::HyperspectralRadiation{FT} = HyperspectralRadiation{FT,DIMS}(DATASET)
end

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-13: move all dimensions here
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Global configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct SPACDimension
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of canopy layers"
    DIM_CANOPY::Int = 20
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Number of wavelength bins for NIR"
    DIM_NIR::Int
    "Number of wavelength bins for PAR"
    DIM_PAR::Int
    "Dimension of root layers"
    DIM_ROOT::Int = 4
    "Dimension of SIF wave length bins"
    DIM_SIF::Int
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int
    "Dimension of soil layers"
    DIM_SOIL::Int = 4
    "Dimension of short wave length bins"
    DIM_WL::Int
    "Dimension of xylem elements"
    DIM_XYLEM::Int = 5
end


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

    # Wavelength information
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT} = WaveLengthSet{FT}(DATASET)

    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to use steady state plant hydraulic system"
    STEADY_STATE_HS::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of canopy layers"
    DIM_CANOPY::Int = 20
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Number of wavelength bins for NIR"
    DIM_NIR::Int = length(WLSET.IΛ_NIR)
    "Number of wavelength bins for PAR"
    DIM_PAR::Int = length(WLSET.IΛ_PAR)
    "Dimension of root layers"
    DIM_ROOT::Int = 4
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = length(WLSET.Λ_SIF)
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = length(WLSET.Λ_SIFE)
    "Dimension of soil layers"
    DIM_SOIL::Int = 4
    "Dimension of short wave length bins"
    DIM_WL::Int = length(WLSET.Λ)
    "Dimension of xylem elements"
    DIM_XYLEM::Int = 5

    # Vectors
    "A matrix of characteristic curves"
    MAT_ρ::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]

    # Embedded structures
    "Hyperspectral absorption features of different leaf components"
    LHA::HyperspectralAbsorption{FT} = HyperspectralAbsorption{FT,DIM_WL}(DATASET)
    "Downwelling shortwave radiation reference spectrum"
    RAD_SW_REF::HyperspectralRadiation{FT} = HyperspectralRadiation{FT,DIM_WL}(DATASET)
end

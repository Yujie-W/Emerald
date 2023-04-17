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
Base.@kwdef struct SPACConfiguration{FT,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}#,SIZE_WL4}
    # Wavelength information
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}

    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to use steady state plant hydraulic system"
    STEADY_STATE_HS::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Vectors
    "A matrix of characteristic curves"
    MAT_ρ::Matrix{FT}#SMatrix{DIM_WL,4,FT,SIZE_WL4}

    # Embedded structures
    "Hyperspectral absorption features of different leaf components"
    LHA::HyperspectralAbsorption{FT,DIM_WL}
    "Downwelling shortwave radiation reference spectrum"
    RAD_SW_REF::HyperspectralRadiation{FT,DIM_WL}
end

SPACConfiguration{FT}(gcf::GeneralConfiguration, dims::SPACDimension) where {FT} = (
    return SPACConfiguration{FT,dims.DIM_NIR,dims.DIM_PAR,dims.DIM_SIF,dims.DIM_SIFE,dims.DIM_WL}(#,dims.DIM_WL*4}(
                WLSET      = WaveLengthSet{FT}(gcf, dims),
                MAT_ρ      = FT[read_nc(gcf.DATASET, "GSV_1") read_nc(gcf.DATASET, "GSV_2") read_nc(gcf.DATASET, "GSV_3") read_nc(gcf.DATASET, "GSV_4")],
                LHA        = HyperspectralAbsorption{FT,dims.DIM_WL}(gcf.DATASET),
                RAD_SW_REF = HyperspectralRadiation{FT,dims.DIM_WL}(gcf.DATASET),
    )
);

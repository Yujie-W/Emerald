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
Base.@kwdef struct SPACConfiguration{FT,DIMS}
    # Wavelength information
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT}

    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to use steady state plant hydraulic system"
    STEADY_STATE_HS::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Vectors
    "A matrix of characteristic curves"
    MAT_ρ::Matrix{FT}

    # Embedded structures
    "Hyperspectral absorption features of different leaf components"
    LHA::HyperspectralAbsorption{FT,DIMS}
    "Downwelling shortwave radiation reference spectrum"
    RAD_SW_REF::HyperspectralRadiation{FT,DIMS}
end

SPACConfiguration{FT,DIMS}(gcf::GeneralConfiguration) where {FT,DIMS} = (
    return SPACConfiguration{FT,DIMS}(
                WLSET      = WaveLengthSet{FT}(gcf,DIMS),
                MAT_ρ      = FT[read_nc(gcf.DATASET, "GSV_1") read_nc(gcf.DATASET, "GSV_2") read_nc(gcf.DATASET, "GSV_3") read_nc(gcf.DATASET, "GSV_4")],
                LHA        = HyperspectralAbsorption{FT,DIMS}(gcf.DATASET),
                RAD_SW_REF = HyperspectralRadiation{FT,DIMS}(gcf.DATASET),
    )
);

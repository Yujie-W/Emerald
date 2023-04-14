#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#     2022-Apr-14: add DIM_WL to HyperspectralAbsorption type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralAbsorption{FT,DIM_WL}
    # Constant features
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::Vector{FT}
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::Vector{FT}
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT}
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT}
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT}
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT}
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT}
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT}
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT}
    "Specific absorption coefficients of PS I and II `[-]`"
    K_PS::Vector{FT}
    "Refractive index `[-]`"
    NR::Vector{FT}
end

HyperspectralAbsorption{FT,DIM_WL}(dset::String) where {FT,DIM_WL} = (
    return HyperspectralAbsorption{FT,DIM_WL}(
                K_ANT   = read_nc(dset, "K_ANT"),
                K_BROWN = read_nc(dset, "K_BROWN"),
                K_CAB   = read_nc(dset, "K_CAB"),
                K_CAR_V = read_nc(dset, "K_CAR_V"),
                K_CAR_Z = read_nc(dset, "K_CAR_Z"),
                K_CBC   = read_nc(dset, "K_CBC"),
                K_H₂O   = read_nc(dset, "K_H₂O"),
                K_LMA   = read_nc(dset, "K_LMA"),
                K_PRO   = read_nc(dset, "K_PRO"),
                K_PS    = read_nc(dset, "K_PS"),
                NR      = read_nc(dset, "NR")
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-10: refactor the structure with renamed fields
#     2021-Aug-10: add a constructor within the structure to avoid external initialization
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#     2022-Apr-14: add DIM_WL to WaveLengthSet type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct WaveLengthSet{FT,DIM_WL}
    # Constants
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT} = FT[700, 2500]
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT} = FT[400, 750]
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT} = FT[640, 850]
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT} = FT[400, 750]
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT}
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::Vector{FT}
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::Vector{FT}

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int} = findall( WL_NIR[1] .<= Λ .<= WL_NIR[2] )
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int} = findall( WL_PAR[1] .<= Λ .<= WL_PAR[2] )
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int} = findall( WL_SIF[1] .<= Λ .<= WL_SIF[2] )
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int} = findall( WL_SIFE[1] .<= Λ .<= WL_SIFE[2] )

    # Constants based on the ones above
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT} = Λ_UPPER .- Λ_LOWER
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT} = Λ[IΛ_PAR]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT} = Λ[IΛ_SIFE]
end

WaveLengthSet{FT}(dset::String) where {FT} =  (
    _nwl = size_nc(dset, "WL");

    return WaveLengthSet{FT,_nwl}(
                Λ       = read_nc(dset, "WL"),
                Λ_LOWER = read_nc(dset, "WL_LOWER"),
                Λ_UPPER = read_nc(dset, "WL_UPPER"),
    )
);


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

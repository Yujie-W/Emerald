#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-13: add state struct to save SPAC configurations
#     2023-Apr-13: move Φ_PHOTON, RAD_SW_REF from MultiLayerSPAC
#     2023-Apr-13: move APAR_CAR from leaf structs
#     2023-Jun-13: add trace gasses as fields
#     2023-Jun-16: move field WLSET from HyperspectralMLCanopy
#     2023-Jun-16: move all fields DIM_* to SPACConfiguration
#     2023-Jun-20: move LHA from SPAC canopy
#     2023-Jun-20: move fields Θ_AZI, Θ_INCL, Θ_INCL_BNDS, _1_AZI, _COS²_Θ_INCL, _COS_Θ_INCL_AZI, and _COS²_Θ_INCL_AZI from spac canopy
#     2023-Jun-20: move fields α_CLM, α_FITTING, and MAT_ρ from soil albedo
#     2023-Jul-06: add field PRESCRIBE_AIR
#     2023-Aug-27: add field ALLOW_LEAF_CONDENSATION
#     2023-Sep-07: add field ALLOW_LEAF_SHEDDING, ENABLE_SOIL_EVAPORATION, and T_CLM
#     2023-Sep-11: add option update legacy to the SPAC configuration
#     2023-Sep-11: add option KR_THRESHOLD to the SPAC configuration
#     2023-Sep-11: add fields ENABLE_ENERGY_BUDGET, ENABLE_PLANT_HYDRAULICS, and ENABLE_SOIL_WATER_BUDGET
#     2023-Sep-14: make sure the configuration struct is consistent with the DATASET
#     2023-Sep-18: add field Φ_SIF_WL for a new feature
#     2023-Sep-19: add field Φ_SIF_CUTOFF and Φ_SIF_RESCALE for a new feature to force SIF wavelength to be lower than the excitation wavelength
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Global configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACConfiguration{FT}
    # Debug mode
    "Whether to print debug information"
    DEBUG::Bool = false

    # Turn on/off features
    "Allow leaf condensation"
    ALLOW_LEAF_CONDENSATION::Bool = false
    "Allow leaf shedding"
    ALLOW_LEAF_SHEDDING::Bool = false
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Enable drought legacy effect"
    ENABLE_DROUGHT_LEGACY::Bool = false
    "Enable energy balance (t_on)"
    ENABLE_ENERGY_BUDGET::Bool = true
    "Enable plant hydraulics (p_on)"
    ENABLE_PLANT_HYDRAULICS::Bool = true
    "Enable soil water budget (θ_on)"
    ENABLE_SOIL_WATER_BUDGET::Bool = true
    "Enable soil evaporation"
    ENABLE_SOIL_EVAPORATION::Bool = false
    "Whether to acclimate leaf Vcmax and Jmax TD"
    T_CLM::Bool = true
    "Whether to use CLM soil albedo scheme"
    α_CLM::Bool = false
    "Whether to fit the data from broadband to hyperspectral"
    α_FITTING::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true
    "How SIF wavelength cutoff is handled (0 for no cut off, 1 for sharp cut off, and 2 for sigmoid used in SCOPE)"
    Φ_SIF_CUTOFF::Int = 0
    "Rescale SIF fluorescence to wavelength lower than the excitation wavelength"
    Φ_SIF_RESCALE::Bool = true
    "Whether to partition the SIF PDF based the wavelength"
    Φ_SIF_WL::Bool = true

    # Prescribe parameters
    "Prescribe air layer information such as partial pressures"
    PRESCRIBE_AIR::Bool = true

    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # WaveLengthSet
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT} = WaveLengthSet{FT}(DATASET = DATASET)
    "Reference Spetra"
    SPECTRA::ReferenceSpectra{FT} = ReferenceSpectra{FT}(DATASET = DATASET)

    # Dimensions of the spac system
    "Dimension of air layers"
    DIM_AIR::Int = 20
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 12
    "Number of wavelength bins for NIR"
    DIM_NIR::Int = length(WLSET.IΛ_NIR)
    "Dimension of root layers"
    DIM_ROOT::Int = 5
    "Number of wavelength bins for PAR"
    DIM_PAR::Int = length(WLSET.IΛ_PAR)
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = length(WLSET.IΛ_SIF)
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = length(WLSET.IΛ_SIFE)
    "Dimension of soil layers"
    DIM_SOIL::Int = 4
    "Dimension of short wave length bins"
    DIM_WL::Int = length(WLSET.Λ)

    # Constant values used to configurate the thresholds
    "Threshold of the critical pressure or flow that trigger a remainder of conductance"
    KR_THRESHOLD::FT = 0.001

    # Canopy geometry related angles
    "Mean azimuth angles `[°]`"
    Θ_AZI::Vector{FT} = collect(FT, range(0, 360; length=DIM_AZI+1))[1:end-1] .+ 360 / DIM_AZI / 2
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT} = FT[ collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1] collect(FT, range(0, 90; length=DIM_INCL+1))[2:end] ]
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT} = FT[ (Θ_INCL_BNDS[_i,1] + Θ_INCL_BNDS[_i,2]) / 2 for _i in 1:DIM_INCL ]

    # Trace gas information
    "Trace gas air"
    TRACE_AIR::TraceGasAir = TraceGasAir{FT}()
    "Trace gas CH₄"
    TRACE_CH₄::TraceGasCH₄ = TraceGasCH₄{FT}()
    "Trace gas CO₂"
    TRACE_CO₂::TraceGasCO₂ = TraceGasCO₂{FT}()
    "Trace gas H₂O"
    TRACE_H₂O::TraceGasH₂O = TraceGasH₂O{FT}()
    "Trace gas N₂"
    TRACE_N₂::TraceGasN₂ = TraceGasN₂{FT}()
    "Trace gas O₂"
    TRACE_O₂::TraceGasO₂ = TraceGasO₂{FT}()

    # Cache variables
    "Ones with the length of Θ_AZI"
    _1_AZI::Vector{FT} = ones(FT, DIM_AZI)
    "Cosine of Θ_AZI"
    _COS_Θ_AZI::Vector{FT} = cosd.(Θ_AZI)
    "Square of cosine of Θ_INCL"
    _COS²_Θ_INCL::Vector{FT} = cosd.(Θ_INCL) .^ 2
    "Square of cosine of Θ_INCL at different azimuth angles"
    _COS²_Θ_INCL_AZI::Matrix{FT} = (cosd.(Θ_INCL) .^ 2) * _1_AZI'
end

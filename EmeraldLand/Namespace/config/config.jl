# This file contains the global configuration of the SPAC system
# The first step for running any code should be to define a config struct

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-13: add state struct to save SPAC configurations
#     2023-Apr-13: add fields Φ_PHOTON, RAD_SW_REF, and APAR_CAR
#     2023-Jun-13: add trace gasses as fields
#     2023-Jun-16: add fields WLSET, DIM_*
#     2023-Jun-20: add fields LHA, Θ_AZI, Θ_INCL, Θ_INCL_BNDS, _COS_Θ_INCL_AZI, α_CLM, α_FITTING, and MAT_ρ
#     2023-Jul-06: add field PRESCRIBE_AIR
#     2023-Aug-27: add field ALLOW_LEAF_CONDENSATION
#     2023-Sep-07: add fields ALLOW_LEAF_SHEDDING, and T_CLM
#     2023-Sep-11: add fields ENABLE_DROUGHT_LEGACY, KR_THRESHOLD
#     2023-Sep-18: add field Φ_SIF_WL
#     2023-Sep-19: add fields Φ_SIF_CUTOFF, and Φ_SIF_RESCALE
#     2023-Sep-20: add new meta field SPECTRA (WLSET, MAT_ρ, ...)
#     2023-Oct-02: add field MESSAGE_LEVEL
#     2023-Oct-05: add field ENABLE_SIF
#     2023-Oct-14: add field ENABLE_REF
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
    "Message level (0 for no, 1 for progress bar, and 2 for ind)"
    MESSAGE_LEVEL::Int = 0

    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # Reference spectra
    "Reference Spetra"
    SPECTRA::ReferenceSpectra{FT} = ReferenceSpectra{FT}(DATASET = DATASET)

    # features related to the leaf optics and fluorescence
    "Whether to compute canopy reflectance"
    ENABLE_REF::Bool = true
    "Whether to compute fluorescence"
    ENABLE_SIF::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true
    "How SIF wavelength cutoff is handled (0 for no cut off, 1 for sharp cut off, and 2 for sigmoid used in SCOPE)"
    Φ_SIF_CUTOFF::Int = 0
    "Rescale SIF fluorescence to wavelength lower than the excitation wavelength"
    Φ_SIF_RESCALE::Bool = true
    "Whether to partition the SIF PDF based the wavelength"
    Φ_SIF_WL::Bool = true

    # features related to plant hydraulics
    "Dimension of xylem slices of leaf, stem, and root; xylem capaciatance of stem and root"
    DIM_XYLEM::Int = 5
    "Threshold of the critical pressure or flow that trigger a remainder of conductance"
    KR_THRESHOLD::FT = 0.001
    "Whether to run the model at steady state mode"
    STEADY_STATE_FLOW::Bool = true

    # features related to canopy photosynthesis
    # Note
    #     1. to use the hyperspectral mode, set both to true
    #     2. to use the broadband mode with sunlit/shaded fractions, set SUNLIT_FRACTION to true and SUNLIT_ANGLES to false
    #     3. to use the broadband mode with one leaf model, set both to false (ppar_sunlit and ppar_shaded will be set to be the same)
    #     4. to use big leaf model, TODO item
    "Whether to partition the sunlit fraction into different inclination and azimuth angles (if false, use float for sunlit fraction)"
    SUNLIT_ANGLES::Bool = true
    "Whether to partition the canopy into sunlit and shaded fractions (if false, use one leaf model)"
    SUNLIT_FRACTION::Bool = true

    # Settings related to the canopy geometry angles (inclination and azimuth settings)
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Mean azimuth angles `[°]`"
    Θ_AZI::Vector{FT} = collect(FT, range(0, 360; length=DIM_AZI+1))[1:end-1] .+ 360 / DIM_AZI / 2
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT} = FT[ collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1] collect(FT, range(0, 90; length=DIM_INCL+1))[2:end] ]
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT} = FT[ (Θ_INCL_BNDS[i,1] + Θ_INCL_BNDS[i,2]) / 2 for i in 1:DIM_INCL ]






    # Turn on/off features
    "Allow leaf condensation"
    ALLOW_LEAF_CONDENSATION::Bool = false
    "Allow leaf shedding"
    ALLOW_LEAF_SHEDDING::Bool = false
    "Enable drought legacy effect"
    ENABLE_DROUGHT_LEGACY::Bool = false
    "Whether to acclimate leaf Vcmax and Jmax TD"
    T_CLM::Bool = true
    "Whether to use CLM soil albedo scheme"
    α_CLM::Bool = false
    "Whether to fit the data from broadband to hyperspectral"
    α_FITTING::Bool = true

    # Prescribe parameters
    "Prescribe air layer information such as partial pressures"
    PRESCRIBE_AIR::Bool = true

    # Trace gas information
    "Trace gas air"
    TRACE_AIR::TraceGasAir{FT} = TraceGasAir{FT}()
    "Trace gas CH₄"
    TRACE_CH₄::TraceGasCH₄{FT} = TraceGasCH₄{FT}()
    "Trace gas CO₂"
    TRACE_CO₂::TraceGasCO₂{FT} = TraceGasCO₂{FT}()
    "Trace gas H₂O"
    TRACE_H₂O::TraceGasH₂O{FT} = TraceGasH₂O{FT}()
    "Trace gas N₂"
    TRACE_N₂::TraceGasN₂{FT} = TraceGasN₂{FT}()
    "Trace gas O₂"
    TRACE_O₂::TraceGasO₂{FT} = TraceGasO₂{FT}()
end;

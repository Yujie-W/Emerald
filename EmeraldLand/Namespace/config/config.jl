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
#     2024-Feb-27: add field HOT_SPOT and SOIL_ALBEDO
#     2024-Feb-27: add constructor function to create the configuration using customized settings
#     2024-Feb-28: add field VERTICAL_BIO
#     2024-Jul-22: add field ENABLE_CHEMICAL_ENERGY
#     2024-Jul-30: do not bin PPAR if DIM_PPAR_BINS is nothing
#     2024-Jul-31: add rate constants fields (constants for PSI, PSII, and combined)
#     2024-Aug-01: add field ENABLE_KD_TD and FIX_ETA_TD
#     2024-Aug-05: set ENABLE_DROUGHT_LEGACY to true by default (add corresponding functions in PlantHydraulics module)
#     2024-Aug-06: set root disconnection threshold to 0.5 loss of root conductance (to avoid numerical issues)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Global configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACConfiguration{FT}
    #
    # Debug mode
    #
    "Message level (0 for no, 1 for progress bar, and 2 for ind)"
    MESSAGE_LEVEL::Int = 0

    #
    # Dataset (to use the reference spectra and many others)
    #
    "Data name of the JLD2 file"
    DATASET::String = OLD_PHI_2021
    "JLD2 file name"
    JLD2_FILE::String = LAND_ARTIFACT
    "Reference Spetra"
    SPECTRA::ReferenceSpectra{FT} = ReferenceSpectra{FT}(JLD2_FILE, DATASET)

    #
    # Canopy optics
    #
    "Whether to compute canopy reflectance"
    ENABLE_REF::Bool = true
    "Hot spot parameter"
    HOT_SPOT::FT = 0.05
    "Soil albedo method"
    SOIL_ALBEDO::Union{SoilAlbedoPrescribe, SoilAlbedoBroadbandCLM, SoilAlbedoBroadbandCLIMA, SoilAlbedoHyperspectralCLM, SoilAlbedoHyperspectralCLIMA} = SoilAlbedoBroadbandCLIMA()
    "Vertical distribution of leaf biophysical properties (if false, run leaf_spectra! only once)"
    VERTICAL_BIO::Bool = false

    #
    # Canopy structure
    #
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

    #
    # Fluorescence
    #
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

    #
    # Photosynthesis
    #
    "Number of sunlit PPAR bins for all the layers (to speed up the computation; 0 for one leaf model)"
    DIM_PPAR_BINS::Union{Int,Nothing} = nothing
    "Enable the chemical energy related to photosynthesis and respiration"
    ENABLE_CHEMICAL_ENERGY::Bool = true
    "Enable K_D rate constant temperature dependency"
    ENABLE_KD_TD::Bool = true
    "Fix the TD of η for Johnson-Berry model"
    FIX_ETA_TD::Bool = true
    "Rate constants for PSI and PSII combined (most for PSII?)"
    PS_RATE_CONSTANTS::PhotosystemsRateConstants{FT} = PhotosystemsRateConstants{FT}()
    "Rate constants for PSI"
    PSI_RATE_CONSTANTS::PhotosystemIRateConstants{FT} = PhotosystemIRateConstants{FT}()
    "Rate constants for PSII"
    PSII_RATE_CONSTANTS::PhotosystemIIRateConstants{FT} = PhotosystemIIRateConstants{FT}()
    "Whether to acclimate leaf Vcmax and Jmax TD"
    T_CLM::Bool = true
    "Number of temperature memory (hours for hourly simulation)"
    T_CLM_N::Int = 240

    #
    # Plant hydraulics
    #
    "Allow leaf condensation"
    ALLOW_LEAF_CONDENSATION::Bool = false
    "Allow leaf shedding in prescibe LAI mode to avoid numerical issues"
    ALLOW_LEAF_SHEDDING::Bool = true
    "Dimension of xylem slices of leaf, stem, and root; xylem capaciatance of stem and root"
    DIM_XYLEM::Int = 5
    "Enable drought legacy effect"
    ENABLE_DROUGHT_LEGACY::Bool = true
    "Threshold of the critical pressure or flow that trigger root disconnection"
    KR_ROOT_DISCONNECTION::FT = 0.5
    "Threshold of the critical pressure or flow that trigger a remainder of conductance"
    KR_THRESHOLD::FT = 0.001
    "Whether to run the model at steady state mode"
    STEADY_STATE_FLOW::Bool = true

    #
    # Trace gas
    #
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

    #
    # Prescribe parameters
    #
    "Prescribe air layer information such as partial pressures"
    PRESCRIBE_AIR::Bool = true
end;

SPACConfiguration(FT::DataType; dataset::String = OLD_PHI_2021, jld2_file::String = LAND_ARTIFACT, wl_par::Vector = [300,750], wl_par_700::Vector = [300,700]) = SPACConfiguration{FT}(
            JLD2_FILE = jld2_file,
            DATASET   = dataset,
            SPECTRA   = ReferenceSpectra{FT}(jld2_file, dataset; wl_par = wl_par, wl_par_700 = wl_par_700),
);

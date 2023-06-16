#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-13: add state struct to save SPAC configurations
#     2023-Apr-13: move Φ_PHOTON, RAD_SW_REF from MultiLayerSPAC
#     2023-Apr-13: move APAR_CAR from leaf structs
#     2023-Jun-13: add trace gasses as fields
#     2023-Jun-16: move field WLSET from HyperspectralMLCanopy
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Global configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACConfiguration{FT}
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT} = WaveLengthSet{FT}()

    # Dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 12
    "Number of wavelength bins for NIR"
    DIM_NIR::Int = length(WLSET.IΛ_NIR)
    "Number of wavelength bins for PAR"
    DIM_PAR::Int = length(WLSET.IΛ_PAR)
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = length(WLSET.IΛ_SIF)
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = length(WLSET.IΛ_SIFE)
    "Dimension of short wave length bins"
    DIM_WL::Int = length(WLSET.Λ)

    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Embedded structures
    "Downwelling shortwave radiation reference spectrum"
    RAD_SW_REF::HyperspectralRadiation{FT} = HyperspectralRadiation{FT}()

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
end

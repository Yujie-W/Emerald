"""
Method configuration for the SPAC model
"""
Base.@kwdef mutable struct SPACConstants{FT<:AbstractFloat}
    # photosynthesis rate constants
    "Rate constants for PSI and PSII combined (most for PSII?)"
    PS_RATE_CONSTANTS::PhotosystemsRateConstants{FT} = PhotosystemsRateConstants{FT}()
    "Rate constants for PSI"
    PSI_RATE_CONSTANTS::PhotosystemIRateConstants{FT} = PhotosystemIRateConstants{FT}()
    "Rate constants for PSII"
    PSII_RATE_CONSTANTS::PhotosystemIIRateConstants{FT} = PhotosystemIIRateConstants{FT}()

    # spectra
    "Reference Spetra"
    SPECTRA::ReferenceSpectra{FT} = ReferenceSpectra{FT}(LAND_ARTIFACT, OLD_PHI_2021)

    # trace gas
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

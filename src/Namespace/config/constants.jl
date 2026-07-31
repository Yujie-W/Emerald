"""
Method configuration for the SPAC model
"""
Base.@kwdef mutable struct SPACConstants{FT<:AbstractFloat}
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

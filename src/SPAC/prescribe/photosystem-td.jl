"""

    prescribe_ps_td!(config::SPACConfig{FT}; args...) where {FT}

Prescribe the photosystem temperature dependence, given
- `config` Configuration for SPAC, with the arguments:
    - `t_acclim` Moving average temperature to update Vcmax and Jmax temperature dependencies

"""
function prescribe_ps_td!(config::SPACConfig{FT}; args...) where {FT}
    prescribe_ps_td_vcmax!(config, config.METHODS.PS_METHODS.TD_VCMAX_C3; args...);
    prescribe_ps_td_vcmax!(config, config.METHODS.PS_METHODS.TD_VCMAX_C4; args...);
    prescribe_ps_td_jmax!(config, config.METHODS.PS_METHODS.TD_JMAX; args...);

    return nothing
end;

# non-public functions to prescribe Vcmax and Jmax TD
prescribe_ps_td_vcmax!(config::SPACConfig{FT}, vtd::AbstractTemperatureDependency{FT}; args...) where {FT} = nothing;

prescribe_ps_td_vcmax!(config::SPACConfig{FT}, vtd::Union{ArrheniusPeak{FT}, ArrheniusPeak2{FT}, Q10Peak{FT}}; t_acclim::Union{Nothing,Number} = nothing) where {FT} = (
    (; ACCLIMATE_T_VCMAX) = config.FEATURES;

    if ACCLIMATE_T_VCMAX && !isnothing(t_acclim)
        vtd.ΔSV = 668.39 - 1.07 * (t_acclim - T₀(FT));
    end;

    return nothing
);

prescribe_ps_td_jmax!(config::SPACConfig{FT}, jtd::AbstractTemperatureDependency{FT}; args...) where {FT} = nothing;

prescribe_ps_td_jmax!(config::SPACConfig{FT}, jtd::Union{ArrheniusPeak{FT}, ArrheniusPeak2{FT}, Q10Peak{FT}}; t_acclim::Union{Nothing,Number} = nothing) where {FT} = (
    (; ACCLIMATE_T_VCMAX) = config.FEATURES;

    if ACCLIMATE_T_VCMAX && !isnothing(t_acclim)
        jtd.ΔSV = 659.70 - 0.75 * (t_acclim - T₀(FT));
    end;

    return nothing
);

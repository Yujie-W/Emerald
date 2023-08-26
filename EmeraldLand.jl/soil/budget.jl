#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add function for water budget
#     2022-Jun-14: use K_MAX and ΔZ and remove K_REF
#     2022-Jun-14: rescale rain for layer 1
#     2022-Jun-14: use METEO.rain
#     2022-Jun-14: add function for soil energy budget
#     2022-Jun-14: use METEO.rain and METEO.t_precip
#     2022-Jun-14: add net radiation energy to top soil
#     2022-Jun-15: add controller to make sure soil layers do not over saturate
#     2022-Jun-15: merge the soil_water! and soil_energy! to soil_budget!
#     2022-Jun-16: move time stepper controller to SoilPlantAirContinuum.jl
#     2022-Jul-26: fix the unit of rain, mass flow, and root extraction (all in mol s⁻¹)
#     2022-Sep-07: allow soil water oversaturation
#     2023-Mar-27: fix a typo when updating e per layer (should use ΔZ per layer rather than the first layer)
#     2023-Apr-07: fix a typo when updating water content in saturated soil layers
#     2023-Apr-08: make runoff a cumulative value within a time interval
#     2023-Jun-13: add trace gas diffusions
#     2023-Jun-13: add diffusion related water and energy budgets
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Jun-30: use separated function for better readability
#     2023-Jul-06: add DEBUG code block
#
#######################################################################################################################################################################################################
"""
#
    soil_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Update the marginal increase of soil water content and energy per layer, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

#
    soil_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat}

Run soil water and energy budget, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
function soil_budget! end

soil_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    soil_infiltration!(config, spac);
    soil_diffusion!(config, spac);
    soil_source_sink!(config, spac);

    return nothing
);

soil_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    soil_diffusion!(config, spac, δt);
    volume_balance!(config, spac);
    soil_infiltration!(config, spac, δt);
    surface_runoff!(config, spac);

    (; DEBUG) = config;
    (; SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # update soil temperature at each layer (top layer t will be same as _t above)
    for _i in eachindex(LAYERS)
        _slayer = LAYERS[_i];
        _cp_gas = (_slayer.TRACES.n_H₂O * CP_V_MOL(FT) + (_slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂) * CP_D_MOL(FT)) / _slayer.ΔZ;
        _slayer._cp = _slayer.ρ * _slayer.CP + _slayer.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas;
        _slayer.t = _slayer.e / _slayer._cp;

        if DEBUG
            if any(isnan, (_cp_gas, _slayer._cp, _slayer.t))
                @info "Debugging" _cp_gas _slayer._cp _slayer.t;
                error("NaN detected when computing soil temperature at layer $(_i)");
            end;
        end;
    end;

    return nothing
);

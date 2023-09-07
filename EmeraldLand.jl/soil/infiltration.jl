#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-29: copy function out of soil_budget!
#     2023-Jun-30: add DEBUG mode to check for NaNs
#     2023-Jul-06: add info into DEBUG code block
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-07: add integrators for soil water budget
#
#######################################################################################################################################################################################################
"""
#
    soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the infiltration rate between soil layers, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

#
    soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update soil water content and energy per layer, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
function soil_infiltration! end

soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; DEBUG, DIM_SOIL) = config;
    (; METEO, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # update soil k, kd, ψ, and λ_thermal for each soil layer (0.5 for tortuosity factor)
    for _slayer in LAYERS
        _slayer.k          = relative_hydraulic_conductance(_slayer.VC, _slayer.θ) * _slayer.VC.K_MAX * relative_viscosity(_slayer.t) / _slayer.ΔZ;
        _slayer.ψ          = soil_ψ_25(_slayer.VC, _slayer.θ; oversaturation = true) * relative_surface_tension(_slayer.t);
        _slayer._kd        = 0.5 * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / _slayer.ΔZ;
        _slayer._λ_thermal = (_slayer.Λ_THERMAL + _slayer.θ * Λ_THERMAL_H₂O(FT)) / _slayer.ΔZ;
        _slayer.∂e∂t       = 0;
        _slayer.∂n∂t      .= 0;
        _slayer.∂θ∂t       = 0;

        if DEBUG
            if any(isnan, (_slayer.k, _slayer.ψ, _slayer._kd, _slayer._λ_thermal))
                @info "Debugging" _slayer.θ;
                @info "Debugging" _slayer.k _slayer.ψ _slayer._kd _slayer._λ_thermal;
                error("NaN detected in soil_infiltration! when computing soil layer hydraulic and thermal properties");
            end;
        end;
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[1].ΔZ;
    LAYERS[1].∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    LAYERS[1].∂e∂t += SOIL.ALBEDO.r_net_lw + SOIL.ALBEDO.r_net_sw;

    if DEBUG
        if any(isnan, (LAYERS[1].∂θ∂t, LAYERS[1].∂e∂t))
            @info "Debugging" LAYERS[1].∂θ∂t LAYERS[1].∂e∂t;
            error("NaN detected in soil_infiltration! when computing top layer water and energy budget from rain and radiation");
        end;
    end;

    for _i in 1:DIM_SOIL-1
        SOIL._k[_i]         = 1 / (2 / LAYERS[_i].k + 2 / LAYERS[_i+1].k);
        SOIL._δψ[_i]        = LAYERS[_i].ψ - LAYERS[_i+1].ψ + ρg_MPa(FT) * (LAYERS[_i].Z - LAYERS[_i+1].Z);
        SOIL._q[_i]         = SOIL._k[_i] * SOIL._δψ[_i];
        SOIL._λ_thermal[_i] = 1 / (2 / LAYERS[_i]._λ_thermal + 2 / LAYERS[_i+1]._λ_thermal);
        SOIL._δt[_i]        = LAYERS[_i].t - LAYERS[_i+1].t;
        SOIL._q_thermal[_i] = SOIL._λ_thermal[_i] * SOIL._δt[_i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (SOIL._q[_i] > 0) && (LAYERS[_i+1].θ >= LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (SOIL._q[_i] < 0) && (LAYERS[_i].θ >= LAYERS[_i].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if both layers are oversaturated, move the oversaturated part from lower layer to upper layer
        if (LAYERS[_i].θ >= LAYERS[_i].VC.Θ_SAT) && (LAYERS[_i+1].θ > LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = -1 * (LAYERS[_i+1].θ - LAYERS[_i+1].VC.Θ_SAT) * LAYERS[_i+1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);
        end;

        if DEBUG
            if any(isnan, (SOIL._k[_i], SOIL._δψ[_i], SOIL._q[_i], SOIL._λ_thermal[_i], SOIL._δt[_i], SOIL._q_thermal[_i]))
                @info "Debugging" SOIL._k[_i] SOIL._δψ[_i] SOIL._q[_i] SOIL._λ_thermal[_i] SOIL._δt[_i] SOIL._q_thermal[_i];
                error("NaN detected in soil_infiltration! when computing water and energy flow between soil layers");
            end;
        end;

        LAYERS[_i  ].∂θ∂t -= SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂θ∂t += SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q_thermal[_i];
        LAYERS[_i+1].∂e∂t += SOIL._q_thermal[_i];
        LAYERS[_i  ].∂e∂t -= SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
        LAYERS[_i+1].∂e∂t += SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;

        if DEBUG
            if any(isnan, (LAYERS[_i].∂θ∂t, LAYERS[_i+1].∂θ∂t, LAYERS[_i].∂e∂t, LAYERS[_i+1].∂e∂t))
                @info "Debugging" LAYERS[_i].∂θ∂t LAYERS[_i+1].∂θ∂t LAYERS[_i].∂e∂t LAYERS[_i+1].∂e∂t;
                error("NaN detected in soil_infiltration! when computing water and energy budget from flow between soil layers");
            end;
        end;
    end;

    return nothing
);

soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; ALLOW_SOIL_EVAPORATION, DEBUG) = config;
    (; METEO, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # run the water transport (condensation + mass flow)
    for _i in eachindex(LAYERS)
        _slayer = LAYERS[_i];

        # account for evaporation and condensation to/from the air space
        if ALLOW_SOIL_EVAPORATION
            _ps = saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000);
            _δθ_v = (_slayer.TRACES.n_H₂O / _slayer.ΔZ - _ps * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R(FT) * _slayer.t)) * M_H₂O(FT) / ρ_H₂O(FT);

            _slayer.θ += _δθ_v;
            _slayer.e += _δθ_v * ρ_H₂O(FT) * CP_V(FT) * _slayer.t; # this energy is transferred from/to air, so use CP_V
            _slayer.e += _δθ_v * ρ_H₂O(FT) * latent_heat_vapor(_slayer.t);

            if DEBUG
                if any(isnan, (_ps, _δθ_v, _slayer.θ, _slayer.e))
                    @info "Debugging" _ps _δθ_v _slayer.θ _slayer.e;
                    error("NaN detected in soil_infiltration! when computing water budget from condensation and mass flow");
                end;
            end;
        end;

        # account for mass flow
        _slayer.θ += _slayer.∂θ∂t * δt;
        _slayer.e += _slayer.∂e∂t * δt / _slayer.ΔZ;
        _slayer.∫∂w∂t_out -= _slayer.∂θ∂t * _slayer.ΔZ * δt *  ρ_H₂O(FT) / M_H₂O(FT);
        if _i == 1
            _slayer.∫∂w∂t_out += METEO.rain * δt;
        end;

        if DEBUG
            if any(isnan, (_slayer.θ, _slayer.e))
                @info "Debugging" _slayer.θ _slayer.e;
                error("NaN detected in soil_infiltration! when computing water budget from condensation and mass flow");
            end;
        end;
    end;

    return nothing
);

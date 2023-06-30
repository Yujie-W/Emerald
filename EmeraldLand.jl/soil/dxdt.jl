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
#
#######################################################################################################################################################################################################
"""
This function have two major features:
- Compute the marginal change of soil water content and energy
- Update soil water content and energy without over-saturating or draining the soil

"""
function soil_budget! end

"""

    soil_budget!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat}

Update the marginal increase of soil water content and energy per layer, given
- `spac` `MultiLayerSPAC` SPAC
- `config` Configuration for `MultiLayerSPAC`

"""
soil_budget!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat} = (
    (; AIR, METEO, ROOTS, ROOTS_INDEX, SOIL) = spac;
    (; DIM_SOIL, TRACE_AIR, TRACE_CH₄, TRACE_CO₂, TRACE_H₂O, TRACE_N₂, TRACE_O₂) = config;
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
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[1].ΔZ;
    LAYERS[1].∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    LAYERS[1].∂e∂t += SOIL.ALBEDO.r_net_lw + SOIL.ALBEDO.r_net_sw;
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

        LAYERS[_i  ].∂θ∂t -= SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂θ∂t += SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q_thermal[_i];
        LAYERS[_i+1].∂e∂t += SOIL._q_thermal[_i];
        LAYERS[_i  ].∂e∂t -= SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
        LAYERS[_i+1].∂e∂t += SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    # update diffusion rate among layers
    _alayer = AIR[1];
    _factor = 2 * LAYERS[1]._kd;
    _v_gas = max(0, LAYERS[1].VC.Θ_SAT - LAYERS[1].θ);

    if _v_gas > 0
        LAYERS[1].∂n∂t[1] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_CH₄, TRACE_AIR) * (LAYERS[1].TRACES.n_CH₄ / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_CH₄ / (GAS_R(FT) * _alayer.t));
        LAYERS[1].∂n∂t[2] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_CO₂, TRACE_AIR) * (LAYERS[1].TRACES.n_CO₂ / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_CO₂ / (GAS_R(FT) * _alayer.t));
        LAYERS[1].∂n∂t[3] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_H₂O, TRACE_AIR) * (LAYERS[1].TRACES.n_H₂O / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_H₂O / (GAS_R(FT) * _alayer.t));
        LAYERS[1].∂n∂t[4] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_N₂ , TRACE_AIR) * (LAYERS[1].TRACES.n_N₂  / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_N₂  / (GAS_R(FT) * _alayer.t));
        LAYERS[1].∂n∂t[5] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_O₂ , TRACE_AIR) * (LAYERS[1].TRACES.n_O₂  / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_O₂  / (GAS_R(FT) * _alayer.t));
    end;

    for _i in 1:DIM_SOIL-1
        # gas diffusion
        _ratei1 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_CH₄, TRACE_AIR);
        _ratei2 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_CO₂, TRACE_AIR);
        _ratei3 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_H₂O, TRACE_AIR);
        _ratei4 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_N₂ , TRACE_AIR);
        _ratei5 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_O₂ , TRACE_AIR);
        _ratej1 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_CH₄, TRACE_AIR);
        _ratej2 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_CO₂, TRACE_AIR);
        _ratej3 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_H₂O, TRACE_AIR);
        _ratej4 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_N₂ , TRACE_AIR);
        _ratej5 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_O₂ , TRACE_AIR);

        if LAYERS[_i]._kd * LAYERS[_i+1]._kd > 0
            _ratio1 = 2 * LAYERS[_i]._kd * _ratei1 * LAYERS[_i+1]._kd * _ratej1 / (LAYERS[_i]._kd * _ratei1 + LAYERS[_i+1]._kd * _ratej1);
            _ratio2 = 2 * LAYERS[_i]._kd * _ratei2 * LAYERS[_i+1]._kd * _ratej2 / (LAYERS[_i]._kd * _ratei2 + LAYERS[_i+1]._kd * _ratej2);
            _ratio3 = 2 * LAYERS[_i]._kd * _ratei3 * LAYERS[_i+1]._kd * _ratej3 / (LAYERS[_i]._kd * _ratei3 + LAYERS[_i+1]._kd * _ratej3);
            _ratio4 = 2 * LAYERS[_i]._kd * _ratei4 * LAYERS[_i+1]._kd * _ratej4 / (LAYERS[_i]._kd * _ratei4 + LAYERS[_i+1]._kd * _ratej4);
            _ratio5 = 2 * LAYERS[_i]._kd * _ratei5 * LAYERS[_i+1]._kd * _ratej5 / (LAYERS[_i]._kd * _ratei5 + LAYERS[_i+1]._kd * _ratej5);
            _drate1 = _ratio1 * (LAYERS[_i].TRACES.n_CH₄ / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_CH₄ / LAYERS[_i+1].ΔZ);
            _drate2 = _ratio2 * (LAYERS[_i].TRACES.n_CO₂ / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_CO₂ / LAYERS[_i+1].ΔZ);
            _drate3 = _ratio3 * (LAYERS[_i].TRACES.n_H₂O / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_H₂O / LAYERS[_i+1].ΔZ);
            _drate4 = _ratio4 * (LAYERS[_i].TRACES.n_N₂  / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_N₂  / LAYERS[_i+1].ΔZ);
            _drate5 = _ratio5 * (LAYERS[_i].TRACES.n_O₂  / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_O₂  / LAYERS[_i+1].ΔZ);
        else
            _ratio1 = _ratio2 = _ratio3 = _ratio4 = _ratio5 = 0;
            _drate1 = _drate2 = _drate3 = _drate4 = _drate5 = 0;
        end;

        LAYERS[_i  ].∂n∂t[1] -= _drate1;
        LAYERS[_i  ].∂n∂t[2] -= _drate2;
        LAYERS[_i  ].∂n∂t[3] -= _drate3;
        LAYERS[_i  ].∂n∂t[4] -= _drate4;
        LAYERS[_i  ].∂n∂t[5] -= _drate5;
        LAYERS[_i+1].∂n∂t[1] += _drate1;
        LAYERS[_i+1].∂n∂t[2] += _drate2;
        LAYERS[_i+1].∂n∂t[3] += _drate3;
        LAYERS[_i+1].∂n∂t[4] += _drate4;
        LAYERS[_i+1].∂n∂t[5] += _drate5;

        # energy transfer related to gas diffusion
        _δe_gas = ((_drate1 + _drate2 + _drate4 + _drate5) * CP_D_MOL(FT) + _drate3 * CP_V_MOL(FT)) * LAYERS[_i].t;
        LAYERS[_i  ].∂e∂t -= _δe_gas;
        LAYERS[_i+1].∂e∂t += _δe_gas;
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) * M_H₂O(FT) / ρ_H₂O(FT) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    return nothing
);

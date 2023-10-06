#######################################################################################################################################################################################################
#
# Changes to the function
# General
#    2023-Jun-29: copy function out of soil_budget!
#    2023-Jun-29: add code to display debug information
#     2023-Jul-06: add info into DEBUG code block
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-07: fix a typo in the concentration calculations
#     2023-Sep-09: fix a typo in the concentration calculations
#     2023-Sep-11: rename ALLOW_SOIL_EVAPORATION to ENABLE_SOIL_EVAPORATION
#
#######################################################################################################################################################################################################
"""
#
    soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the diffusion rate among soil layers, given
- `config` SPAC configuration
- `spac` SPAC model

#
    soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update diffusion rate among soil layers (and thus water and energy budgets), given
- `config` SPAC configuration
- `spac` SPAC model
- `δt` time step

"""
function soil_diffusion! end

soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; ENABLE_SOIL_EVAPORATION, DEBUG) = config;

    if !ENABLE_SOIL_EVAPORATION
        return nothing
    end;

    (; AIR, SOILS) = spac;
    (; DIM_SOIL, TRACE_AIR, TRACE_CH₄, TRACE_CO₂, TRACE_H₂O, TRACE_N₂, TRACE_O₂) = config;

    # update diffusion rate among layers
    _alayer = AIR[1];
    _factor = 2 * SOILS[1]._kd;
    _v_gas = max(0, SOILS[1].VC.Θ_SAT - SOILS[1].θ) * SOILS[1].ΔZ;
    if _v_gas > 0
        SOILS[1].∂n∂t[1] -= _factor * diffusive_coefficient(SOILS[1].t, TRACE_CH₄, TRACE_AIR) * (SOILS[1].state.ns[1] / _v_gas - _alayer.p_CH₄ / (GAS_R(FT) * _alayer.t));
        SOILS[1].∂n∂t[2] -= _factor * diffusive_coefficient(SOILS[1].t, TRACE_CO₂, TRACE_AIR) * (SOILS[1].state.ns[2] / _v_gas - _alayer.p_CO₂ / (GAS_R(FT) * _alayer.t));
        SOILS[1].∂n∂t[3] -= _factor * diffusive_coefficient(SOILS[1].t, TRACE_H₂O, TRACE_AIR) * (SOILS[1].state.ns[3] / _v_gas - _alayer.p_H₂O / (GAS_R(FT) * _alayer.t));
        SOILS[1].∂n∂t[4] -= _factor * diffusive_coefficient(SOILS[1].t, TRACE_N₂ , TRACE_AIR) * (SOILS[1].state.ns[4]  / _v_gas - _alayer.p_N₂  / (GAS_R(FT) * _alayer.t));
        SOILS[1].∂n∂t[5] -= _factor * diffusive_coefficient(SOILS[1].t, TRACE_O₂ , TRACE_AIR) * (SOILS[1].state.ns[5]  / _v_gas - _alayer.p_O₂  / (GAS_R(FT) * _alayer.t));
    end;

    if DEBUG
        if any(isnan, SOILS[1].∂n∂t) || any(isnan, (_factor, _v_gas))
            @info "Debugging" _factor _v_gas SOILS[1].∂n∂t[1] SOILS[1].∂n∂t[2] SOILS[1].∂n∂t[3] SOILS[1].∂n∂t[4] SOILS[1].∂n∂t[5];
            error("NaN in soil_diffusion! at layer 1");
        end;
    end;

    # TODO: did I forget this term? (SOILS[1].ΔZ * _v_gas)
    for i in 1:DIM_SOIL-1
        _v_gas_i = max(0, SOILS[i].VC.Θ_SAT - SOILS[i].θ) * SOILS[i].ΔZ;
        _v_gas_j = max(0, SOILS[i+1].VC.Θ_SAT - SOILS[i+1].θ) * SOILS[i+1].ΔZ;

        # gas diffusion
        _ratei1 = diffusive_coefficient(SOILS[i  ].t, TRACE_CH₄, TRACE_AIR);
        _ratei2 = diffusive_coefficient(SOILS[i  ].t, TRACE_CO₂, TRACE_AIR);
        _ratei3 = diffusive_coefficient(SOILS[i  ].t, TRACE_H₂O, TRACE_AIR);
        _ratei4 = diffusive_coefficient(SOILS[i  ].t, TRACE_N₂ , TRACE_AIR);
        _ratei5 = diffusive_coefficient(SOILS[i  ].t, TRACE_O₂ , TRACE_AIR);
        _ratej1 = diffusive_coefficient(SOILS[i+1].t, TRACE_CH₄, TRACE_AIR);
        _ratej2 = diffusive_coefficient(SOILS[i+1].t, TRACE_CO₂, TRACE_AIR);
        _ratej3 = diffusive_coefficient(SOILS[i+1].t, TRACE_H₂O, TRACE_AIR);
        _ratej4 = diffusive_coefficient(SOILS[i+1].t, TRACE_N₂ , TRACE_AIR);
        _ratej5 = diffusive_coefficient(SOILS[i+1].t, TRACE_O₂ , TRACE_AIR);

        if DEBUG
            if any(isnan, (_v_gas_i, _v_gas_j, _ratei1, _ratei2, _ratei3, _ratei4, _ratei5, _ratej1, _ratej2, _ratej3, _ratej4, _ratej5))
                @info "Debugging" _v_gas_i _v_gas_j _ratei1 _ratei2 _ratei3 _ratei4 _ratei5 _ratej1 _ratej2 _ratej3 _ratej4 _ratej5;
                error("NaN in soil_diffusion! at layers $(i) and $(i+1) when computing diffusion coefficient");
            end;
        end;

        if SOILS[i]._kd * SOILS[i+1]._kd > 0
            _ratio1 = 2 * SOILS[i]._kd * _ratei1 * SOILS[i+1]._kd * _ratej1 / (SOILS[i]._kd * _ratei1 + SOILS[i+1]._kd * _ratej1);
            _ratio2 = 2 * SOILS[i]._kd * _ratei2 * SOILS[i+1]._kd * _ratej2 / (SOILS[i]._kd * _ratei2 + SOILS[i+1]._kd * _ratej2);
            _ratio3 = 2 * SOILS[i]._kd * _ratei3 * SOILS[i+1]._kd * _ratej3 / (SOILS[i]._kd * _ratei3 + SOILS[i+1]._kd * _ratej3);
            _ratio4 = 2 * SOILS[i]._kd * _ratei4 * SOILS[i+1]._kd * _ratej4 / (SOILS[i]._kd * _ratei4 + SOILS[i+1]._kd * _ratej4);
            _ratio5 = 2 * SOILS[i]._kd * _ratei5 * SOILS[i+1]._kd * _ratej5 / (SOILS[i]._kd * _ratei5 + SOILS[i+1]._kd * _ratej5);
            _drate1 = _ratio1 * (SOILS[i].state.ns[1] / _v_gas_i - SOILS[i+1].state.ns[1] / _v_gas_j);
            _drate2 = _ratio2 * (SOILS[i].state.ns[2] / _v_gas_i - SOILS[i+1].state.ns[2] / _v_gas_j);
            _drate3 = _ratio3 * (SOILS[i].state.ns[3] / _v_gas_i - SOILS[i+1].state.ns[3] / _v_gas_j);
            _drate4 = _ratio4 * (SOILS[i].state.ns[4]  / _v_gas_i - SOILS[i+1].state.ns[4]  / _v_gas_j);
            _drate5 = _ratio5 * (SOILS[i].state.ns[5]  / _v_gas_i - SOILS[i+1].state.ns[5]  / _v_gas_j);
        else
            _ratio1 = _ratio2 = _ratio3 = _ratio4 = _ratio5 = 0;
            _drate1 = _drate2 = _drate3 = _drate4 = _drate5 = 0;
        end;

        if DEBUG
            if any(isnan, (_ratio1, _ratio2, _ratio3, _ratio4, _ratio5, _drate1, _drate2, _drate3, _drate4, _drate5))
                @info "Debugging" SOILS[i]._kd SOILS[i+1]._kd _v_gas_i _v_gas_j _ratio1 _ratio2 _ratio3 _ratio4 _ratio5 _drate1 _drate2 _drate3 _drate4 _drate5;
                error("NaN in soil_diffusion! at layers $(i) and $(i+1) when computing diffusion rate");
            end;
        end;

        SOILS[i  ].∂n∂t[1] -= _drate1;
        SOILS[i  ].∂n∂t[2] -= _drate2;
        SOILS[i  ].∂n∂t[3] -= _drate3;
        SOILS[i  ].∂n∂t[4] -= _drate4;
        SOILS[i  ].∂n∂t[5] -= _drate5;
        SOILS[i+1].∂n∂t[1] += _drate1;
        SOILS[i+1].∂n∂t[2] += _drate2;
        SOILS[i+1].∂n∂t[3] += _drate3;
        SOILS[i+1].∂n∂t[4] += _drate4;
        SOILS[i+1].∂n∂t[5] += _drate5;

        if DEBUG
            if any(isnan, (SOILS[i].∂n∂t[1], SOILS[i].∂n∂t[2], SOILS[i].∂n∂t[3], SOILS[i].∂n∂t[4], SOILS[i].∂n∂t[5],
                           SOILS[i+1].∂n∂t[1], SOILS[i+1].∂n∂t[2], SOILS[i+1].∂n∂t[3], SOILS[i+1].∂n∂t[4], SOILS[i+1].∂n∂t[5]))
                @info "Debugging" SOILS[i].∂n∂t[1] SOILS[i].∂n∂t[2] SOILS[i].∂n∂t[3] SOILS[i].∂n∂t[4] SOILS[i].∂n∂t[5];
                @info "Debugging" SOILS[i+1].∂n∂t[1], SOILS[i+1].∂n∂t[2], SOILS[i+1].∂n∂t[3], SOILS[i+1].∂n∂t[4], SOILS[i+1].∂n∂t[5];
                error("NaN in soil_diffusion! at layers $(i) and $(i+1) when computing ∂n∂t");
            end;
        end;

        # energy transfer related to gas diffusion
        _δe_gas = ((_drate1 + _drate2 + _drate4 + _drate5) * CP_D_MOL(FT) + _drate3 * CP_V_MOL(FT)) * SOILS[i].t;
        SOILS[i  ].∂e∂t -= _δe_gas;
        SOILS[i+1].∂e∂t += _δe_gas;

        if DEBUG
            if any(isnan, (_δe_gas, SOILS[i].∂e∂t, SOILS[i+1].∂e∂t))
                @info "Debugging" _δe_gas SOILS[i].∂e∂t SOILS[i+1].∂e∂t;
                error("NaN in soil_diffusion! at layers $(i) and $(i+1) when computing ∂e∂t");
            end;
        end;
    end;

    return nothing
);

soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; ENABLE_SOIL_EVAPORATION, DEBUG) = config;

    if !ENABLE_SOIL_EVAPORATION
        return nothing
    end;

    (; SOILS) = spac;

    # run the diffusion
    for soil in SOILS
        _δθ = max(0, soil.VC.Θ_SAT - soil.θ);
        if _δθ == 0
            soil.state.ns[1] = 0;
            soil.state.ns[2] = 0;
            soil.state.ns[3] = 0;
            soil.state.ns[4]  = 0;
            soil.state.ns[5]  = 0;
        else
            soil.state.ns[1] += soil.∂n∂t[1] * δt;
            soil.state.ns[2] += soil.∂n∂t[2] * δt;
            soil.state.ns[3] += soil.∂n∂t[3] * δt;
            soil.state.ns[4]  += soil.∂n∂t[4] * δt;
            soil.state.ns[5]  += soil.∂n∂t[5] * δt;
        end;

        if DEBUG
            if any(isnan, (_δθ, soil.state.ns[1], soil.state.ns[2], soil.state.ns[3], soil.state.ns[4], soil.state.ns[5]))
                @info "Debugging" _δθ soil.state.ns[1] soil.state.ns[2] soil.state.ns[3] soil.state.ns[4] soil.state.ns[5];
                error("NaN detected in soil_diffusion! dxdt");
            end;
        end;
    end;

    return nothing
);

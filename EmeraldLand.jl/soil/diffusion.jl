#######################################################################################################################################################################################################
#
# Changes to the function
# General
#    2023-Jun-29: copy function out of soil_budget!
#    2023-Jun-29: add code to display debug information
#
#######################################################################################################################################################################################################
"""
#
    soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Compute the diffusion rate among soil layers, given
- `config` SPAC configuration
- `spac` SPAC model

#
    soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat}

Update diffusion rate among soil layers (and thus water and energy budgets), given
- `config` SPAC configuration
- `spac` SPAC model
- `δt` time step

"""
function soil_diffusion! end

soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; DEBUG) = config;
    (; AIR, SOIL) = spac;
    (; DIM_SOIL, TRACE_AIR, TRACE_CH₄, TRACE_CO₂, TRACE_H₂O, TRACE_N₂, TRACE_O₂) = config;
    LAYERS = SOIL.LAYERS;

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

    if DEBUG
        if any(isnan, LAYERS[1].∂n∂t) || any(isnan, (_factor, _v_gas))
            @error "NaN in soil_diffusion! at layer 1";
        end;
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

        if DEBUG
            if any(isnan, (_ratei1, _ratei2, _ratei3, _ratei4, _ratei5, _ratej1, _ratej2, _ratej3, _ratej4, _ratej5))
                @error "NaN in soil_diffusion! at layers $(_i) and $(_i+1) when computing diffusion coefficient";
            end;
        end;

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

        if DEBUG
            if any(isnan, (_ratio1, _ratio2, _ratio3, _ratio4, _ratio5, _drate1, _drate2, _drate3, _drate4, _drate5))
                @error "NaN in soil_diffusion! at layers $(_i) and $(_i+1) when computing diffusion rate";
            end;
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

        if DEBUG
            if any(isnan, (LAYERS[_i].∂n∂t[1], LAYERS[_i].∂n∂t[2], LAYERS[_i].∂n∂t[3], LAYERS[_i].∂n∂t[4], LAYERS[_i].∂n∂t[5],
                           LAYERS[_i+1].∂n∂t[1], LAYERS[_i+1].∂n∂t[2], LAYERS[_i+1].∂n∂t[3], LAYERS[_i+1].∂n∂t[4], LAYERS[_i+1].∂n∂t[5]))
                @error "NaN in soil_diffusion! at layers $(_i) and $(_i+1) when computing ∂n∂t";
            end;
        end;

        # energy transfer related to gas diffusion
        _δe_gas = ((_drate1 + _drate2 + _drate4 + _drate5) * CP_D_MOL(FT) + _drate3 * CP_V_MOL(FT)) * LAYERS[_i].t;
        LAYERS[_i  ].∂e∂t -= _δe_gas;
        LAYERS[_i+1].∂e∂t += _δe_gas;

        if DEBUG
            if any(isnan, (_δe_gas, LAYERS[_i].∂e∂t, LAYERS[_i+1].∂e∂t))
                @error "NaN in soil_diffusion! at layers $(_i) and $(_i+1) when computing ∂e∂t";
            end;
        end;
    end;

    return nothing
);

soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    (; DEBUG) = config;
    (; SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # run the diffusion
    for _slayer in LAYERS
        _δθ = max(0, _slayer.VC.Θ_SAT - _slayer.θ);
        if _δθ == 0
            _slayer.TRACES.n_CH₄ = 0;
            _slayer.TRACES.n_CO₂ = 0;
            _slayer.TRACES.n_H₂O = 0;
            _slayer.TRACES.n_N₂  = 0;
            _slayer.TRACES.n_O₂  = 0;
        else
            _slayer.TRACES.n_CH₄ += _slayer.∂n∂t[1] * δt;
            _slayer.TRACES.n_CO₂ += _slayer.∂n∂t[2] * δt;
            _slayer.TRACES.n_H₂O += _slayer.∂n∂t[3] * δt;
            _slayer.TRACES.n_N₂  += _slayer.∂n∂t[4] * δt;
            _slayer.TRACES.n_O₂  += _slayer.∂n∂t[5] * δt;
        end;

        if DEBUG
            if any(isnan, (_δθ, _slayer.TRACES.n_CH₄, _slayer.TRACES.n_CO₂, _slayer.TRACES.n_H₂O, _slayer.TRACES.n_N₂, _slayer.TRACES.n_O₂))
                @error "NaN detected in soil_diffusion! dxdt";
            end;
        end;
    end;

    return nothing
);

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#    2023-Jun-30: move function out of soil_budget!
#
#######################################################################################################################################################################################################
"""

    volume_balance!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Balance the air volume in the soil so that pressure is in equilibrium, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

"""
function volume_balance!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}
    (; AIR, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # compute air volume change using ideal gas law (total energy change accordingly)
    println("\nJudge point 9");
    _alayer = AIR[1];
    for _i in length(LAYERS):-1:1
        _slayer = LAYERS[_i];
        _n_air = (_alayer.P_AIR - saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000)) * _slayer.ΔZ * _slayer.θ / (GAS_R(FT) * _slayer.t);
        _n_dry = _slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂;
        _n_rat = (_n_air - _n_dry) / _n_dry;

        @show _n_air _n_dry;

        # change air moles anyway because top soil is exposed to the atmosphere
        if (_i == 1) && (_n_air - _n_dry != 0)
            if _n_dry == 0
                _m_dry = _alayer.n_CH₄ + _alayer.n_CO₂ + _alayer.n_N₂ + _alayer.n_O₂;
                _m_rat = (_n_air - _n_dry) / _m_dry;
                _slayer.TRACES.n_CH₄ += _m_rat * _alayer.n_CH₄;
                _slayer.TRACES.n_CO₂ += _m_rat * _alayer.n_CO₂;
                _slayer.TRACES.n_N₂  += _m_rat * _alayer.n_N₂;
                _slayer.TRACES.n_O₂  += _m_rat * _alayer.n_O₂;
                _slayer.e += (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
                _alayer.n_CH₄ -= _m_rat * _alayer.n_CH₄;
                _alayer.n_CO₂ -= _m_rat * _alayer.n_CO₂;
                _alayer.n_N₂  -= _m_rat * _alayer.n_N₂;
                _alayer.n_O₂  -= _m_rat * _alayer.n_O₂;
                _alayer.e -= (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;

                @info "i == 1, n_dry == 0, n_air - n_dry != 0";
            else
                _slayer.TRACES.n_CH₄ += _n_rat * _slayer.TRACES.n_CH₄;
                _slayer.TRACES.n_CO₂ += _n_rat * _slayer.TRACES.n_CO₂;
                _slayer.TRACES.n_N₂  += _n_rat * _slayer.TRACES.n_N₂;
                _slayer.TRACES.n_O₂  += _n_rat * _slayer.TRACES.n_O₂;
                _slayer.e += (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
                _alayer.n_CH₄ -= _n_rat * _slayer.TRACES.n_CH₄;
                _alayer.n_CO₂ -= _n_rat * _slayer.TRACES.n_CO₂;
                _alayer.n_N₂  -= _n_rat * _slayer.TRACES.n_N₂;
                _alayer.n_O₂  -= _n_rat * _slayer.TRACES.n_O₂;
                _alayer.e -= (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;

                @info "i == 1, n_dry != 0, n_air - n_dry != 0";
            end;

            @show _slayer.TRACES.n_CH₄ _slayer.TRACES.n_CO₂ _slayer.TRACES.n_H₂O _slayer.TRACES.n_N₂ _slayer.TRACES.n_O₂;

        elseif _n_air - _n_dry != 0
            _tlayer = LAYERS[_i-1];
            _δvt = _tlayer.ΔZ * max(0, _tlayer.VC.Θ_SAT - _tlayer.θ);

            @show _δvt;

            if _n_dry == 0
                if _δvt > 0
                    _m_dry = _tlayer.TRACES.n_CH₄ + _tlayer.TRACES.n_CO₂ + _tlayer.TRACES.n_N₂ + _tlayer.TRACES.n_O₂;
                    _m_rat = (_n_air - _n_dry) / _m_dry;
                    _slayer.TRACES.n_CH₄ += _m_rat * _tlayer.TRACES.n_CH₄;
                    _slayer.TRACES.n_CO₂ += _m_rat * _tlayer.TRACES.n_CO₂;
                    _slayer.TRACES.n_N₂  += _m_rat * _tlayer.TRACES.n_N₂;
                    _slayer.TRACES.n_O₂  += _m_rat * _tlayer.TRACES.n_O₂;
                    _slayer.e += (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
                    _tlayer.TRACES.n_CH₄ -= _m_rat * _tlayer.TRACES.n_CH₄;
                    _tlayer.TRACES.n_CO₂ -= _m_rat * _tlayer.TRACES.n_CO₂;
                    _tlayer.TRACES.n_N₂  -= _m_rat * _tlayer.TRACES.n_N₂;
                    _tlayer.TRACES.n_O₂  -= _m_rat * _tlayer.TRACES.n_O₂;
                    _tlayer.e -= (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;

                    @info "i != 1, n_dry == 0, n_air - n_dry != 0, δvt > 0";
                else
                    _θ_fc_up = soil_θ(_tlayer.VC, -1 * ρ_H₂O(FT) * GRAVITY(FT) * _tlayer.ΔZ / relative_surface_tension(_tlayer.t));
                    _p_dry = _alayer.P_AIR - saturation_vapor_pressure(_tlayer.t);
                    _max_n_t = _p_dry * (_tlayer.VC.Θ_SAT - _θ_fc_up) * _tlayer.ΔZ / (GAS_R(FT) * _tlayer.t);

                    # suck water directly from upper layer
                    @info "Drawing water from upper layer 0";
                    _v_up = min(_max_n_t, abs(_n_air - _n_dry)) * GAS_R(FT) * _tlayer.t / _p_dry;

                    @show _max_n_t _p_dry _θ_fc_up _v_up;

                    _slayer.θ += _v_up / _slayer.ΔZ;
                    _tlayer.θ -= _v_up / _tlayer.ΔZ;
                    _slayer.e += _v_up * ρ_H₂O(FT) * CP_L(FT) * _slayer.t;
                    _tlayer.e -= _v_up * ρ_H₂O(FT) * CP_L(FT) * _slayer.t;

                    @show _slayer.θ _tlayer.θ _slayer.e _tlayer.e;

                    # if the lower air is releasing air to upper layer, replace the air with water
                    if _n_air - _n_dry < 0
                        @info "Releasing air to upper layer 0";
                        _m_rat = min(_max_n_t, abs(_n_air - _n_dry)) / _n_dry;

                        @show _m_rat;

                        _slayer.TRACES.n_CH₄ -= _m_rat * _slayer.TRACES.n_CH₄;
                        _slayer.TRACES.n_CO₂ -= _m_rat * _slayer.TRACES.n_CO₂;
                        _slayer.TRACES.n_N₂  -= _m_rat * _slayer.TRACES.n_N₂;
                        _slayer.TRACES.n_O₂  -= _m_rat * _slayer.TRACES.n_O₂;
                        _slayer.e -= _p_dry * _v_up / (GAS_R(FT) * _slayer.t) * CP_D_MOL(FT) * _slayer.t;
                        _tlayer.TRACES.n_CH₄ += _m_rat * _slayer.TRACES.n_CH₄;
                        _tlayer.TRACES.n_CO₂ += _m_rat * _slayer.TRACES.n_CO₂;
                        _tlayer.TRACES.n_N₂  += _m_rat * _slayer.TRACES.n_N₂;
                        _tlayer.TRACES.n_O₂  += _m_rat * _slayer.TRACES.n_O₂;
                        _tlayer.e += _p_dry * _v_up / (GAS_R(FT) * _slayer.t) * CP_D_MOL(FT) * _slayer.t;
                    end;
                end;
            else
                # change air moles among layers if the up layer is not saturated
                if _δvt > 0
                    _slayer.TRACES.n_CH₄ += _n_rat * _slayer.TRACES.n_CH₄;
                    _slayer.TRACES.n_CO₂ += _n_rat * _slayer.TRACES.n_CO₂;
                    _slayer.TRACES.n_N₂  += _n_rat * _slayer.TRACES.n_N₂;
                    _slayer.TRACES.n_O₂  += _n_rat * _slayer.TRACES.n_O₂;
                    _slayer.e += (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
                    _tlayer.TRACES.n_CH₄ -= _n_rat * _slayer.TRACES.n_CH₄;
                    _tlayer.TRACES.n_CO₂ -= _n_rat * _slayer.TRACES.n_CO₂;
                    _tlayer.TRACES.n_N₂  -= _n_rat * _slayer.TRACES.n_N₂;
                    _tlayer.TRACES.n_O₂  -= _n_rat * _slayer.TRACES.n_O₂;
                    _tlayer.e -= (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
                # replace air moles with water if the up layer is not saturated
                else
                    _θ_fc_up = soil_θ(_tlayer.VC, -1 * ρ_H₂O(FT) * GRAVITY(FT) * _tlayer.ΔZ / relative_surface_tension(_tlayer.t));
                    _p_dry = _alayer.P_AIR - saturation_vapor_pressure(_tlayer.t);
                    _max_n_t = _p_dry * (_tlayer.VC.Θ_SAT - _θ_fc_up) * _tlayer.ΔZ / (GAS_R(FT) * _tlayer.t);

                    # suck water directly from upper layer
                    @info "Drawing water from upper layer";
                    _v_up = min(_max_n_t, abs(_n_air - _n_dry)) * GAS_R(FT) * _tlayer.t / _p_dry;

                    @show _max_n_t _p_dry _θ_fc_up _v_up;

                    _slayer.θ += _v_up / _slayer.ΔZ;
                    _tlayer.θ -= _v_up / _tlayer.ΔZ;
                    _slayer.e += _v_up * ρ_H₂O(FT) * CP_L(FT) * _slayer.t;
                    _tlayer.e -= _v_up * ρ_H₂O(FT) * CP_L(FT) * _slayer.t;

                    @show _slayer.θ _tlayer.θ _slayer.e _tlayer.e;

                    # if the lower air is releasing air to upper layer, replace the air with water
                    if _n_air - _n_dry < 0
                        @info "Releasing air to upper layer";
                        _m_rat = min(_max_n_t, abs(_n_air - _n_dry)) / _n_dry;

                        @show _m_rat;

                        _slayer.TRACES.n_CH₄ -= _m_rat * _slayer.TRACES.n_CH₄;
                        _slayer.TRACES.n_CO₂ -= _m_rat * _slayer.TRACES.n_CO₂;
                        _slayer.TRACES.n_N₂  -= _m_rat * _slayer.TRACES.n_N₂;
                        _slayer.TRACES.n_O₂  -= _m_rat * _slayer.TRACES.n_O₂;
                        _slayer.e -= _p_dry * _v_up / (GAS_R(FT) * _slayer.t) * CP_D_MOL(FT) * _slayer.t;
                        _tlayer.TRACES.n_CH₄ += _m_rat * _slayer.TRACES.n_CH₄;
                        _tlayer.TRACES.n_CO₂ += _m_rat * _slayer.TRACES.n_CO₂;
                        _tlayer.TRACES.n_N₂  += _m_rat * _slayer.TRACES.n_N₂;
                        _tlayer.TRACES.n_O₂  += _m_rat * _slayer.TRACES.n_O₂;
                        _tlayer.e += _p_dry * _v_up / (GAS_R(FT) * _slayer.t) * CP_D_MOL(FT) * _slayer.t;
                    end;
                end;
            end;

            @show _slayer.TRACES.n_CH₄ _slayer.TRACES.n_CO₂ _slayer.TRACES.n_H₂O _slayer.TRACES.n_N₂ _slayer.TRACES.n_O₂;


        end;
    end;

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#    2023-Jun-30: move function out of soil_budget!
#
#######################################################################################################################################################################################################
"""

    surface_runoff!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute surface runoff, given
- `config` spac configuration
- `spac` spac model

"""
function surface_runoff!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; DEBUG) = config;
    (; SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # compute surface runoff
    if LAYERS[1].θ > LAYERS[1].VC.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp = LAYERS[1].CP * LAYERS[1].ρ + LAYERS[1].θ * ρ_H₂O(FT) * CP_L(FT);
        _t  = LAYERS[1].e / _cp;
        _runoff = (LAYERS[1].θ - LAYERS[1].VC.Θ_SAT) * LAYERS[1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);

        LAYERS[1].θ = LAYERS[1].VC.Θ_SAT;
        LAYERS[1].e -= _runoff / LAYERS[1].ΔZ * CP_L_MOL(FT) * _t;
        SOIL.runoff += _runoff;

        if DEBUG
            if any(isnan, [_cp, _t, _runoff, LAYERS[1].θ, LAYERS[1].e, SOIL.runoff])
                @error "NaN detected in surface_runoff!";
            end;
        end;
    end;

    return nothing
end

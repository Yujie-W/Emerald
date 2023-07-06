#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-30: move function out of soil_budget!
#     2023-Jul-06: sort the order of gas volume balance and water volume balance
#     2023-Jul-06: add DEBUG code block
#     2023-Jul-06: add PRESCRIBE_AIR mode to avoid the errors due to mass balance in air
#
#######################################################################################################################################################################################################
"""

    volume_balance!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Balance the air volume in the soil so that pressure is in equilibrium, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

"""
function volume_balance!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}
    (; DEBUG, PRESCRIBE_AIR) = config;
    (; AIR, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # balance the air volume among soil layers from lower to upper layers
    _alayer = AIR[1];
    for _i in length(LAYERS)-1:-1:1
        # upper layer is _ilayer and lower is _jlayer
        _ilayer = LAYERS[_i];
        _jlayer = LAYERS[_i+1];

        # compute the air moles in the lower layer
        _i_dry = _ilayer.TRACES.n_CH₄ + _ilayer.TRACES.n_CO₂ + _ilayer.TRACES.n_N₂ + _ilayer.TRACES.n_O₂;
        _j_dry = _jlayer.TRACES.n_CH₄ + _jlayer.TRACES.n_CO₂ + _jlayer.TRACES.n_N₂ + _jlayer.TRACES.n_O₂;
        _j_max = (_alayer.P_AIR - saturation_vapor_pressure(_jlayer.t, _jlayer.ψ * 1000000)) * _jlayer.ΔZ * (_jlayer.VC.Θ_SAT - _jlayer.θ) / (GAS_R(FT) * _jlayer.t);

        if DEBUG
            if any(isnan, (_i_dry, _j_dry, _j_max))
                @info "Debugging" _i_dry _j_dry _j_max;
                @error "NaN detected in volume_balance! between layer $(_i) and layer $(_i+1)";
            end;
        end;

        # if _j_dry == _j_max, no air needs to be transferred from/to the upper layer

        # if _j_max > _j_dry and _i_dry > 0, air needs to be transferred from the upper layer
        if (_j_max > _j_dry) && (_i_dry > 0)
            _n_mass = min(_j_max - _j_dry, _i_dry);
            _jlayer.TRACES.n_CH₄ += _n_mass * _ilayer.TRACES.n_CH₄ / _i_dry;
            _jlayer.TRACES.n_CO₂ += _n_mass * _ilayer.TRACES.n_CO₂ / _i_dry;
            _jlayer.TRACES.n_N₂  += _n_mass * _ilayer.TRACES.n_N₂  / _i_dry;
            _jlayer.TRACES.n_O₂  += _n_mass * _ilayer.TRACES.n_O₂  / _i_dry;
            _jlayer.e += _n_mass * CP_D_MOL(FT) * _ilayer.t / _jlayer.ΔZ;
            _ilayer.TRACES.n_CH₄ -= _n_mass * _ilayer.TRACES.n_CH₄ / _i_dry;
            _ilayer.TRACES.n_CO₂ -= _n_mass * _ilayer.TRACES.n_CO₂ / _i_dry;
            _ilayer.TRACES.n_N₂  -= _n_mass * _ilayer.TRACES.n_N₂  / _i_dry;
            _ilayer.TRACES.n_O₂  -= _n_mass * _ilayer.TRACES.n_O₂  / _i_dry;
            _ilayer.e -= _n_mass * CP_D_MOL(FT) * _ilayer.t / _ilayer.ΔZ;

            if DEBUG
                if any(isnan, (_jlayer.TRACES.n_CH₄, _jlayer.TRACES.n_CO₂, _jlayer.TRACES.n_N₂, _jlayer.TRACES.n_O₂, _jlayer.e,
                               _ilayer.TRACES.n_CH₄, _ilayer.TRACES.n_CO₂, _ilayer.TRACES.n_N₂, _ilayer.TRACES.n_O₂, _ilayer.e))
                    @info "Debugging" _jlayer.TRACES.n_CH₄ _jlayer.TRACES.n_CO₂ _jlayer.TRACES.n_N₂ _jlayer.TRACES.n_O₂ _jlayer.e;
                    @info "Debugging" _ilayer.TRACES.n_CH₄ _ilayer.TRACES.n_CO₂ _ilayer.TRACES.n_N₂ _ilayer.TRACES.n_O₂ _ilayer.e;
                    @error "NaN detected in volume_balance! between layer $(_i) and layer $(_i+1) when moving air from upper layer to lower layer";
                end;
            end;
        end;

        # if _j_max > _j_dry but _i_dry == 0, the lower layer will need to suck some water from upper layer to balance the air volume
        # compute the equilibrate mole of air from upper layer when it reaches its field capacity
        if (_j_max > _j_dry) && (_j_dry > 0) && (_i_dry == 0)
            _θ_fc_up = soil_θ(_ilayer.VC, -1 * ρ_H₂O(FT) * GRAVITY(FT) * _ilayer.ΔZ / relative_surface_tension(_ilayer.t));
            _v_mass = min((_j_max - _j_dry) * GAS_R(FT) * _jlayer.t / _alayer.P_AIR, (_ilayer.VC.Θ_SAT - _θ_fc_up) * _ilayer.ΔZ);
            _jlayer.θ += _v_mass / _jlayer.ΔZ;
            _ilayer.θ -= _v_mass / _ilayer.ΔZ;
            _jlayer.e += _v_mass * ρ_H₂O(FT) * CP_L(FT) * _ilayer.t / _jlayer.ΔZ;
            _ilayer.e -= _v_mass * ρ_H₂O(FT) * CP_L(FT) * _ilayer.t / _ilayer.ΔZ;

            if DEBUG
                if any(isnan, (_jlayer.θ, _ilayer.θ, _jlayer.e, _ilayer.e))
                    @info "Debugging" _jlayer.θ _ilayer.θ _jlayer.e _ilayer.e;
                    @error "NaN detected in volume_balance! between layer $(_i) and layer $(_i+1) when moving water from upper layer to lower layer";
                end;
            end;
        end;

        # if _j_max < _j_dry and _j_dry > 0, air needs to be transferred to the upper layer (does not matter whether upper layer is saturated or not)
        if (_j_max < _j_dry) && (_j_dry > 0)
            _n_mass = _j_dry - _j_max;
            _jlayer.TRACES.n_CH₄ -= _n_mass * _jlayer.TRACES.n_CH₄ / _j_dry;
            _jlayer.TRACES.n_CO₂ -= _n_mass * _jlayer.TRACES.n_CO₂ / _j_dry;
            _jlayer.TRACES.n_N₂  -= _n_mass * _jlayer.TRACES.n_N₂  / _j_dry;
            _jlayer.TRACES.n_O₂  -= _n_mass * _jlayer.TRACES.n_O₂  / _j_dry;
            _jlayer.e -= _n_mass * CP_D_MOL(FT) * _jlayer.t / _jlayer.ΔZ;
            _ilayer.TRACES.n_CH₄ += _n_mass * _jlayer.TRACES.n_CH₄ / _j_dry;
            _ilayer.TRACES.n_CO₂ += _n_mass * _jlayer.TRACES.n_CO₂ / _j_dry;
            _ilayer.TRACES.n_N₂  += _n_mass * _jlayer.TRACES.n_N₂  / _j_dry;
            _ilayer.TRACES.n_O₂  += _n_mass * _jlayer.TRACES.n_O₂  / _j_dry;
            _ilayer.e += _n_mass * CP_D_MOL(FT) * _jlayer.t / _ilayer.ΔZ;

            if DEBUG
                if any(isnan, (_jlayer.TRACES.n_CH₄, _jlayer.TRACES.n_CO₂, _jlayer.TRACES.n_N₂, _jlayer.TRACES.n_O₂, _jlayer.e,
                               _ilayer.TRACES.n_CH₄, _ilayer.TRACES.n_CO₂, _ilayer.TRACES.n_N₂, _ilayer.TRACES.n_O₂, _ilayer.e))
                    @info "Debugging" _jlayer.TRACES.n_CH₄ _jlayer.TRACES.n_CO₂ _jlayer.TRACES.n_N₂ _jlayer.TRACES.n_O₂ _jlayer.e;
                    @info "Debugging" _ilayer.TRACES.n_CH₄ _ilayer.TRACES.n_CO₂ _ilayer.TRACES.n_N₂ _ilayer.TRACES.n_O₂ _ilayer.e;
                    @error "NaN detected in volume_balance! between layer $(_i) and layer $(_i+1) when moving air from lower layer to upper layer";
                end;
            end;
        end;
    end;

    # balance the air volume between top soil and atmosphere
    _slayer = LAYERS[1];
    _s_dry = _slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂;
    _a_dry = _alayer.n_CH₄ + _alayer.n_CO₂ + _alayer.n_N₂ + _alayer.n_O₂;
    _s_max = (_alayer.P_AIR - saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000)) * _slayer.ΔZ * (_slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R(FT) * _slayer.t);

    if DEBUG
        if any(isnan, (_s_dry, _a_dry, _s_max))
            @info "Debugging" _s_dry _a_dry _s_max;
            @error "NaN detected in volume_balance! when moving air from top soil from/to atmosphere";
        end;
    end;

    # if soil air is not saturated, it can absorb more air from the atmosphere
    if _s_max > _s_dry
        _n_mass = min(_s_max - _s_dry, _a_dry);
        _slayer.TRACES.n_CH₄ += _n_mass * _alayer.n_CH₄ / _a_dry;
        _slayer.TRACES.n_CO₂ += _n_mass * _alayer.n_CO₂ / _a_dry;
        _slayer.TRACES.n_N₂  += _n_mass * _alayer.n_N₂  / _a_dry;
        _slayer.TRACES.n_O₂  += _n_mass * _alayer.n_O₂  / _a_dry;
        _slayer.e += _n_mass * CP_D_MOL(FT) * _alayer.t / _slayer.ΔZ;
        if !PRESCRIBE_AIR
            _alayer.n_CH₄ -= _n_mass * _alayer.n_CH₄ / _a_dry;
            _alayer.n_CO₂ -= _n_mass * _alayer.n_CO₂ / _a_dry;
            _alayer.n_N₂  -= _n_mass * _alayer.n_N₂  / _a_dry;
            _alayer.n_O₂  -= _n_mass * _alayer.n_O₂  / _a_dry;
            _alayer.e -= _n_mass * CP_D_MOL(FT) * _alayer.t;
        end;

        if DEBUG
            if any(isnan, (_slayer.TRACES.n_CH₄, _slayer.TRACES.n_CO₂, _slayer.TRACES.n_N₂, _slayer.TRACES.n_O₂, _slayer.e,
                           _alayer.n_CH₄, _alayer.n_CO₂, _alayer.n_N₂, _alayer.n_O₂, _alayer.e))
                @info "Debugging" _n_mass _a_dry;
                @info "Debugging" _slayer.TRACES.n_CH₄ _slayer.TRACES.n_CO₂ _slayer.TRACES.n_N₂ _slayer.TRACES.n_O₂ _slayer.e;
                @info "Debugging" _alayer.n_CH₄ _alayer.n_CO₂ _alayer.n_N₂ _alayer.n_O₂ _alayer.e;
                @error "NaN detected in volume_balance! when moving air from atmosphere to top soil";
            end;
        end;
    elseif _s_max < _s_dry
        _n_mass = _s_dry - _s_max;
        _slayer.TRACES.n_CH₄ -= _n_mass * _slayer.TRACES.n_CH₄ / _s_dry;
        _slayer.TRACES.n_CO₂ -= _n_mass * _slayer.TRACES.n_CO₂ / _s_dry;
        _slayer.TRACES.n_N₂  -= _n_mass * _slayer.TRACES.n_N₂  / _s_dry;
        _slayer.TRACES.n_O₂  -= _n_mass * _slayer.TRACES.n_O₂  / _s_dry;
        _slayer.e -= _n_mass * CP_D_MOL(FT) * _slayer.t / _slayer.ΔZ;
        if !PRESCRIBE_AIR
            _alayer.n_CH₄ += _n_mass * _slayer.TRACES.n_CH₄ / _s_dry;
            _alayer.n_CO₂ += _n_mass * _slayer.TRACES.n_CO₂ / _s_dry;
            _alayer.n_N₂  += _n_mass * _slayer.TRACES.n_N₂  / _s_dry;
            _alayer.n_O₂  += _n_mass * _slayer.TRACES.n_O₂  / _s_dry;
            _alayer.e += _n_mass * CP_D_MOL(FT) * _slayer.t;
        end;

        if DEBUG
            if any(isnan, (_slayer.TRACES.n_CH₄, _slayer.TRACES.n_CO₂, _slayer.TRACES.n_N₂, _slayer.TRACES.n_O₂, _slayer.e,
                           _alayer.n_CH₄, _alayer.n_CO₂, _alayer.n_N₂, _alayer.n_O₂, _alayer.e))
                @info "Debugging" _slayer.TRACES.n_CH₄ _slayer.TRACES.n_CO₂ _slayer.TRACES.n_N₂ _slayer.TRACES.n_O₂ _slayer.e;
                @info "Debugging" _alayer.n_CH₄ _alayer.n_CO₂ _alayer.n_N₂ _alayer.n_O₂ _alayer.e;
                @error "NaN detected in volume_balance! when moving air from top soil to atmosphere";
            end;
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
        _cp_gas = (LAYERS[1].TRACES.n_H₂O * CP_V_MOL(FT) + (LAYERS[1].TRACES.n_CH₄ + LAYERS[1].TRACES.n_CO₂ + LAYERS[1].TRACES.n_N₂ + LAYERS[1].TRACES.n_O₂) * CP_D_MOL(FT)) / LAYERS[1].ΔZ;
        _cp = LAYERS[1].CP * LAYERS[1].ρ + LAYERS[1].θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas;
        _t  = LAYERS[1].e / _cp;
        _runoff = (LAYERS[1].θ - LAYERS[1].VC.Θ_SAT) * LAYERS[1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);

        LAYERS[1].θ = LAYERS[1].VC.Θ_SAT;
        LAYERS[1].e -= _runoff / LAYERS[1].ΔZ * CP_L_MOL(FT) * _t;
        SOIL.runoff += _runoff;

        if DEBUG
            if any(isnan, (_cp, _t, _runoff, LAYERS[1].θ, LAYERS[1].e, SOIL.runoff))
                @info "Debugging" _cp _t _runoff LAYERS[1].θ LAYERS[1].e SOIL.runoff;
                @error "NaN detected in surface_runoff!";
            end;
        end;
    end;

    return nothing
end

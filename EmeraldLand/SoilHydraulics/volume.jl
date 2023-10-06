#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-30: move function out of soil_budget!
#     2023-Jul-06: sort the order of gas volume balance and water volume balance
#     2023-Jul-06: add DEBUG code block
#     2023-Jul-06: add PRESCRIBE_AIR mode to avoid the errors due to mass balance in air
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-11: rename ALLOW_SOIL_EVAPORATION to ENABLE_SOIL_EVAPORATION
#
#######################################################################################################################################################################################################
"""

    volume_balance!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Balance the air volume in the soil so that pressure is in equilibrium, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

"""
function volume_balance!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; ENABLE_SOIL_EVAPORATION, DEBUG, PRESCRIBE_AIR) = config;

    if !ENABLE_SOIL_EVAPORATION
        return nothing
    end;

    (; AIR, SOILS) = spac;

    # balance the air volume among soil layers from lower to upper layers
    _alayer = AIR[1];
    for i in length(SOILS)-1:-1:1
        # upper layer is _ilayer and lower is _jlayer
        _ilayer = SOILS[i];
        _jlayer = SOILS[i+1];

        # compute the air moles in the lower layer
        _i_dry = _ilayer.state.ns[1] + _ilayer.state.ns[2] + _ilayer.state.ns[4] + _ilayer.state.ns[5];
        _j_dry = _jlayer.state.ns[1] + _jlayer.state.ns[2] + _jlayer.state.ns[4] + _jlayer.state.ns[5];
        _j_max = (_alayer.P_AIR - saturation_vapor_pressure(_jlayer.t, _jlayer.ψ * 1000000)) * _jlayer.ΔZ * (_jlayer.VC.Θ_SAT - _jlayer.θ) / (GAS_R(FT) * _jlayer.t);

        if DEBUG
            if any(isnan, (_i_dry, _j_dry, _j_max))
                @info "Debugging" _i_dry _j_dry _j_max;
                error("NaN detected in volume_balance! between layer $(i) and layer $(i+1)");
            end;
        end;

        # if _j_dry == _j_max, no air needs to be transferred from/to the upper layer

        # if _j_max > _j_dry and _i_dry > 0, air needs to be transferred from the upper layer
        if (_j_max > _j_dry) && (_i_dry > 0)
            _n_mass = min(_j_max - _j_dry, _i_dry);
            _jlayer.state.ns[1] += _n_mass * _ilayer.state.ns[1] / _i_dry;
            _jlayer.state.ns[2] += _n_mass * _ilayer.state.ns[2] / _i_dry;
            _jlayer.state.ns[4]  += _n_mass * _ilayer.state.ns[4]  / _i_dry;
            _jlayer.state.ns[5]  += _n_mass * _ilayer.state.ns[5]  / _i_dry;
            _jlayer.Σe += _n_mass * CP_D_MOL(FT) * _ilayer.t / _jlayer.ΔZ;
            _ilayer.state.ns[1] -= _n_mass * _ilayer.state.ns[1] / _i_dry;
            _ilayer.state.ns[2] -= _n_mass * _ilayer.state.ns[2] / _i_dry;
            _ilayer.state.ns[4]  -= _n_mass * _ilayer.state.ns[4]  / _i_dry;
            _ilayer.state.ns[5]  -= _n_mass * _ilayer.state.ns[5]  / _i_dry;
            _ilayer.Σe -= _n_mass * CP_D_MOL(FT) * _ilayer.t / _ilayer.ΔZ;

            if DEBUG
                if any(isnan, (_jlayer.state.ns[1], _jlayer.state.ns[2], _jlayer.state.ns[4], _jlayer.state.ns[5], _jlayer.Σe,
                               _ilayer.state.ns[1], _ilayer.state.ns[2], _ilayer.state.ns[4], _ilayer.state.ns[5], _ilayer.Σe))
                    @info "Debugging" _jlayer.state.ns[1] _jlayer.state.ns[2] _jlayer.state.ns[4] _jlayer.state.ns[5] _jlayer.Σe;
                    @info "Debugging" _ilayer.state.ns[1] _ilayer.state.ns[2] _ilayer.state.ns[4] _ilayer.state.ns[5] _ilayer.Σe;
                    error("NaN detected in volume_balance! between layer $(i) and layer $(i+1) when moving air from upper layer to lower layer");
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
            _jlayer.Σe += _v_mass * ρ_H₂O(FT) * CP_L(FT) * _ilayer.t / _jlayer.ΔZ;
            _ilayer.Σe -= _v_mass * ρ_H₂O(FT) * CP_L(FT) * _ilayer.t / _ilayer.ΔZ;

            if DEBUG
                if any(isnan, (_jlayer.θ, _ilayer.θ, _jlayer.Σe, _ilayer.Σe))
                    @info "Debugging" _jlayer.θ _ilayer.θ _jlayer.Σe _ilayer.Σe;
                    error("NaN detected in volume_balance! between layer $(i) and layer $(i+1) when moving water from upper layer to lower layer");
                end;
            end;
        end;

        # if _j_max < _j_dry and _j_dry > 0, air needs to be transferred to the upper layer (does not matter whether upper layer is saturated or not)
        if (_j_max < _j_dry) && (_j_dry > 0)
            _n_mass = _j_dry - _j_max;
            _jlayer.state.ns[1] -= _n_mass * _jlayer.state.ns[1] / _j_dry;
            _jlayer.state.ns[2] -= _n_mass * _jlayer.state.ns[2] / _j_dry;
            _jlayer.state.ns[4]  -= _n_mass * _jlayer.state.ns[4]  / _j_dry;
            _jlayer.state.ns[5]  -= _n_mass * _jlayer.state.ns[5]  / _j_dry;
            _jlayer.Σe -= _n_mass * CP_D_MOL(FT) * _jlayer.t / _jlayer.ΔZ;
            _ilayer.state.ns[1] += _n_mass * _jlayer.state.ns[1] / _j_dry;
            _ilayer.state.ns[2] += _n_mass * _jlayer.state.ns[2] / _j_dry;
            _ilayer.state.ns[4]  += _n_mass * _jlayer.state.ns[4]  / _j_dry;
            _ilayer.state.ns[5]  += _n_mass * _jlayer.state.ns[5]  / _j_dry;
            _ilayer.Σe += _n_mass * CP_D_MOL(FT) * _jlayer.t / _ilayer.ΔZ;

            if DEBUG
                if any(isnan, (_jlayer.state.ns[1], _jlayer.state.ns[2], _jlayer.state.ns[4], _jlayer.state.ns[5], _jlayer.Σe,
                               _ilayer.state.ns[1], _ilayer.state.ns[2], _ilayer.state.ns[4], _ilayer.state.ns[5], _ilayer.Σe))
                    @info "Debugging" _jlayer.state.ns[1] _jlayer.state.ns[2] _jlayer.state.ns[4] _jlayer.state.ns[5] _jlayer.Σe;
                    @info "Debugging" _ilayer.state.ns[1] _ilayer.state.ns[2] _ilayer.state.ns[4] _ilayer.state.ns[5] _ilayer.Σe;
                    error("NaN detected in volume_balance! between layer $(i) and layer $(i+1) when moving air from lower layer to upper layer");
                end;
            end;
        end;
    end;

    # balance the air volume between top soil and atmosphere
    soil = SOILS[1];
    _s_dry = soil.state.ns[1] + soil.state.ns[2] + soil.state.ns[4] + soil.state.ns[5];
    _a_dry = _alayer.n_CH₄ + _alayer.n_CO₂ + _alayer.n_N₂ + _alayer.n_O₂;
    _s_max = (_alayer.P_AIR - saturation_vapor_pressure(soil.t, soil.auxil.ψ * 1000000)) * soil.auxil.δz * (soil.VC.Θ_SAT - soil.θ) / (GAS_R(FT) * soil.t);

    if DEBUG
        if any(isnan, (_s_dry, _a_dry, _s_max))
            @info "Debugging" _s_dry _a_dry _s_max;
            error("NaN detected in volume_balance! when moving air from top soil from/to atmosphere");
        end;
    end;

    # if soil air is not saturated, it can absorb more air from the atmosphere
    if _s_max > _s_dry
        _n_mass = min(_s_max - _s_dry, _a_dry);
        soil.state.ns[1] += _n_mass * _alayer.n_CH₄ / _a_dry;
        soil.state.ns[2] += _n_mass * _alayer.n_CO₂ / _a_dry;
        soil.state.ns[4]  += _n_mass * _alayer.n_N₂  / _a_dry;
        soil.state.ns[5]  += _n_mass * _alayer.n_O₂  / _a_dry;
        soil.state.Σe += _n_mass * CP_D_MOL(FT) * _alayer.t / soil.auxil.δz;
        if !PRESCRIBE_AIR
            _alayer.n_CH₄ -= _n_mass * _alayer.n_CH₄ / _a_dry;
            _alayer.n_CO₂ -= _n_mass * _alayer.n_CO₂ / _a_dry;
            _alayer.n_N₂  -= _n_mass * _alayer.n_N₂  / _a_dry;
            _alayer.n_O₂  -= _n_mass * _alayer.n_O₂  / _a_dry;
            _alayer.Σe -= _n_mass * CP_D_MOL(FT) * _alayer.t;
        end;

        if DEBUG
            if any(isnan, (soil.state.ns[1], soil.state.ns[2], soil.state.ns[4], soil.state.ns[5], soil.state.Σe,
                           _alayer.n_CH₄, _alayer.n_CO₂, _alayer.n_N₂, _alayer.n_O₂, _alayer.Σe))
                @info "Debugging" _n_mass _a_dry;
                @info "Debugging" soil.state.ns[1] soil.state.ns[2] soil.state.ns[4] soil.state.ns[5] soil.state.Σe;
                @info "Debugging" _alayer.n_CH₄ _alayer.n_CO₂ _alayer.n_N₂ _alayer.n_O₂ _alayer.Σe;
                error("NaN detected in volume_balance! when moving air from atmosphere to top soil");
            end;
        end;
    elseif _s_max < _s_dry
        _n_mass = _s_dry - _s_max;
        soil.state.ns[1] -= _n_mass * soil.state.ns[1] / _s_dry;
        soil.state.ns[2] -= _n_mass * soil.state.ns[2] / _s_dry;
        soil.state.ns[4]  -= _n_mass * soil.state.ns[4]  / _s_dry;
        soil.state.ns[5]  -= _n_mass * soil.state.ns[5]  / _s_dry;
        soil.state.Σe -= _n_mass * CP_D_MOL(FT) * soil.t / soil.auxil.δz;
        if !PRESCRIBE_AIR
            _alayer.n_CH₄ += _n_mass * soil.state.ns[1] / _s_dry;
            _alayer.n_CO₂ += _n_mass * soil.state.ns[2] / _s_dry;
            _alayer.n_N₂  += _n_mass * soil.state.ns[4]  / _s_dry;
            _alayer.n_O₂  += _n_mass * soil.state.ns[5]  / _s_dry;
            _alayer.Σe += _n_mass * CP_D_MOL(FT) * soil.t;
        end;

        if DEBUG
            if any(isnan, (soil.state.ns[1], soil.state.ns[2], soil.state.ns[4], soil.state.ns[5], soil.state.Σe,
                           _alayer.n_CH₄, _alayer.n_CO₂, _alayer.n_N₂, _alayer.n_O₂, _alayer.Σe))
                @info "Debugging" soil.state.ns[1] soil.state.ns[2] soil.state.ns[4] soil.state.ns[5] soil.state.Σe;
                @info "Debugging" _alayer.n_CH₄ _alayer.n_CO₂ _alayer.n_N₂ _alayer.n_O₂ _alayer.Σe;
                error("NaN detected in volume_balance! when moving air from top soil to atmosphere");
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
#     2023-Sep-07: add integrators for soil water budget
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
    (; SOIL_BULK, SOILS) = spac;

    # compute surface runoff
    if SOILS[1].state.θ > SOILS[1].state.vc.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp_gas = (SOILS[1].state.ns[3] * CP_V_MOL(FT) + (SOILS[1].state.ns[1] + SOILS[1].state.ns[2] + SOILS[1].state.ns[4] + SOILS[1].state.ns[5]) * CP_D_MOL(FT)) / SOILS[1].auxil.δz;
        _cp = SOILS[1].state.cp * SOILS[1].state.ρ + SOILS[1].state.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas;
        _t  = SOILS[1].state.Σe / _cp;
        _runoff = (SOILS[1].state.θ - SOILS[1].state.vc.Θ_SAT) * SOILS[1].auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);

        SOILS[1].state.θ = SOILS[1].state.vc.Θ_SAT;
        SOILS[1].state.Σe -= _runoff / SOILS[1].auxil.δz * CP_L_MOL(FT) * _t;
        SOIL_BULK.auxil.runoff += _runoff;
    end;

    return nothing
end

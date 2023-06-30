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

    soil_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat}

Run soil water and energy budget, given
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
soil_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    soil_diffusion!(config, spac, δt);

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

    # run the water transport (condensation + mass flow)
    for _slayer in LAYERS
        # account for evaporation and condensation to/from the air space
        _ps = saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000);
        _δθ_v = (_slayer.TRACES.n_H₂O / _slayer.ΔZ - _ps * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R(FT) * _slayer.t)) * M_H₂O(FT) / ρ_H₂O(FT);

        @show _ps _δθ_v;

        _slayer.θ += _δθ_v;
        _slayer.e += _δθ_v * ρ_H₂O(FT) * CP_L(FT) * _slayer.t;
        _slayer.e += _δθ_v * ρ_H₂O(FT) * latent_heat_vapor(_slayer.t);

        # account for mass flow
        _slayer.θ += _slayer.∂θ∂t * δt;
        _slayer.e += _slayer.∂e∂t * δt / _slayer.ΔZ;
    end;

    # compute surface runoff
    if LAYERS[1].θ > LAYERS[1].VC.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp = LAYERS[1].CP * LAYERS[1].ρ + LAYERS[1].θ * ρ_H₂O(FT) * CP_L(FT);
        _t  = LAYERS[1].e / _cp;
        _runoff = (LAYERS[1].θ - LAYERS[1].VC.Θ_SAT) * LAYERS[1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);

        @show _cp _t _runoff;

        LAYERS[1].θ = LAYERS[1].VC.Θ_SAT;
        LAYERS[1].e -= _runoff / LAYERS[1].ΔZ * CP_L_MOL(FT) * _t;
        SOIL.runoff += _runoff;
    end;

    # update soil temperature at each layer (top layer t will be same as _t above)
    for _slayer in LAYERS
        _cp_gas = (_slayer.TRACES.n_H₂O * CP_V_MOL(FT) + (_slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂) * CP_D_MOL(FT)) / _slayer.ΔZ;
        _slayer._cp = _slayer.ρ * _slayer.CP + _slayer.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas;
        _slayer.t = _slayer.e / _slayer._cp;

        @show _cp_gas _slayer._cp _slayer.t;
    end;

    return nothing
);

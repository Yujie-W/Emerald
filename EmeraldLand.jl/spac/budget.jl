#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-15: add function to make sure the soil layers do not over saturate or drain
#     2022-Jun-18: move function from SoilHydraulics.jl to SoilPlantAirContinuum.jl
#     2022-Jun-18: add controller for soil and leaf temperatures
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Aug-31: add controller for leaf stomatal conductance
#     2022-Sep-07: remove soil oversaturation controller, and add a Δθ <= 0.01 controller
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2023-Mar-27: add controller for trunk and branch temperatures
#     2023-Jun-13: add soil gas energy into soil e when computing combined cp
#     2023-Jun-15: add judge for root connection
#     2023-Aug-27: add DEBUG controller
#
#######################################################################################################################################################################################################
"""

    function adjusted_time(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}, δt::FT; t_on::Bool = true, θ_on::Bool = true) where {FT<:AbstractFloat}

Return adjusted time that soil does not over saturate or drain, given
- `spac` `MultiLayerSPAC` SPAC
- `config` Configuration for `MultiLayerSPAC`
- `δt` Time step
- `t_on` If true, plant energy budget is on (set false to run sensitivity analysis or prescribing mode)
- `θ_on` If true, soil water budget is on (set false to run sensitivity analysis or prescribing mode)

"""
function adjusted_time(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}, δt::FT; t_on::Bool = true, θ_on::Bool = true) where {FT<:AbstractFloat}
    (; DEBUG) = config;
    (; BRANCHES, LEAVES, SOIL, TRUNK) = spac;

    _δt_1 = δt;

    # make sure each layer does not drain (allow for oversaturation), and θ change is less than 0.01
    _δt_2 = _δt_1;
    if θ_on
        for _slayer in SOIL.LAYERS
            _δt_2 = min(FT(0.01) / abs(_slayer.∂θ∂t), _δt_2);
            if _slayer.∂θ∂t < 0
                _δt_dra = (_slayer.VC.Θ_RES - _slayer.θ) / _slayer.∂θ∂t;
                _δt_2 = min(_δt_dra, _δt_2);
            end;

            if DEBUG
                if any(isnan, (_δt_2, _slayer.∂θ∂t, _slayer.θ))
                    @info "Debugging" _δt_2 _slayer.∂θ∂t _slayer.θ;
                    error("NaN in adjusted_time at layer 1");
                end;
            end;
        end;
    end;

    # make sure soil temperatures do not change more than 1 K per time step
    _δt_3 = _δt_2;
    if t_on
        for _slayer in SOIL.LAYERS
            _cp_gas = (_slayer.TRACES.n_H₂O * CP_V_MOL(FT) + (_slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂) * CP_D_MOL(FT)) / _slayer.ΔZ;
            _∂T∂t = _slayer.∂e∂t / (_slayer.ρ * _slayer.CP + _slayer.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas);
            _δt_3 = min(1 / abs(_∂T∂t), _δt_3);

            if DEBUG
                if any(isnan, (_δt_3, _∂T∂t, _cp_gas))
                    @info "Debugging" _δt_3 _∂T∂t _cp_gas;
                    error("NaN in adjusted_time at layer 2");
                end;
            end;
        end;
    end;

    # make sure trunk temperatures do not change more than 1 K per time step
    _δt_4 = _δt_3;
    if t_on && spac._root_connection
        _∂T∂t = TRUNK.∂e∂t / (CP_L_MOL(FT) * sum(TRUNK.HS.v_storage));
        _δt_4 = min(1 / abs(_∂T∂t), _δt_4);

        if DEBUG
            if any(isnan, (_δt_4, _∂T∂t))
                @info "Debugging" _δt_4 _∂T∂t;
                error("NaN in adjusted_time at trunk");
            end;
        end;
    end;

    # make sure branch stem temperatures do not change more than 1 K per time step
    _δt_5 = _δt_4;
    if t_on && spac._root_connection
        for _branch in BRANCHES
            _∂T∂t = _branch.∂e∂t / (CP_L_MOL(FT) * sum(_branch.HS.v_storage));
            _δt_5 = min(1 / abs(_∂T∂t), _δt_5);

            if DEBUG
                if any(isnan, (_δt_5, _∂T∂t))
                    @info "Debugging" _δt_5 _∂T∂t;
                    error("NaN in adjusted_time at branch");
                end;
            end;
        end;
    end;

    # make sure leaf temperatures do not change more than 1 K per time step
    _δt_6 = _δt_5;
    if t_on && spac._root_connection
        for _clayer in LEAVES
            _∂T∂t = _clayer.∂e∂t / (_clayer.CP * _clayer.BIO.lma * 10 + CP_L_MOL(FT) * _clayer.HS.v_storage);
            _δt_6 = min(1 / abs(_∂T∂t), _δt_6);

            if DEBUG
                if any(isnan, (_δt_6, _∂T∂t))
                    @info "Debugging" _δt_6 _∂T∂t;
                    error("NaN in adjusted_time at leaf");
                end;
            end;
        end;
    end;

    # make sure leaf stomatal conductances do not change more than 0.06 mol m⁻² s⁻¹
    _δt_7 = _δt_6;
    if spac._root_connection
        for _clayer in LEAVES
            for _∂g∂t in _clayer.∂g∂t_sunlit
                _δt_7 = min(FT(0.06) / abs(_∂g∂t), _δt_7);
            end;
            _δt_7 = min(FT(0.06) / abs(_clayer.∂g∂t_shaded), _δt_7);

            if DEBUG
                if any(isnan, (_δt_7, _clayer.∂g∂t_shaded, mean(_clayer.∂g∂t_sunlit)))
                    @info "Debugging" _δt_7 _clayer.∂g∂t_shaded mean(_clayer.∂g∂t_sunlit);
                    error("NaN in adjusted_time at leaf stomatal conductance");
                end;
            end;
        end;
    end;

    # make sure adjusted time is not nan
    @assert !isnan(_δt_7) "NaN in adjusted_time";

    return (_δt_7, _δt_6, _δt_5, _δt_4, _δt_3, _δt_2, _δt_1)
end


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-18: move function from SoilHydraulics.jl to SoilPlantAirContinuum.jl
#     2022-Jun-18: make it a separate function
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Sep-07: add method to solve for steady state solution
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2022-Nov-18: add option p_on to enable/disable plant flow and pressure profiles
#     2023-Apr-13: add config to function call to steady state function
#     2023-Apr-13: sw and lw radiation moved to METEO
#     2023-Jun-13: add config to parameter list
#     2023-Jun-15: add judge for root connection
#
#######################################################################################################################################################################################################
"""

    time_stepper!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat}
    time_stepper!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}; update::Bool = false) where {FT<:AbstractFloat}

Move forward in time for SPAC with time stepper controller, given
- `spac` `MultiLayerSPAC` SPAC
- `config` Configuration for `MultiLayerSPAC`
- `δt` Time step (if not given, solve for steady state solution)
- `p_on` If true, plant hydraulic flow and pressure profiles will be updated
- `t_on` If true, plant energy budget is on (set false to run sensitivity analysis or prescribing mode)
- `update` If true, update leaf xylem legacy effect
- `θ_on` If true, soil water budget is on (set false to run sensitivity analysis or prescribing mode)

"""
function time_stepper! end

time_stepper!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES, METEO, SOIL) = spac;

    # run the update function until time elapses
    _count = 0;
    _t_res = FT(δt);
    while true
        _count += 1;
        _δts = adjusted_time(spac, config, _t_res; θ_on = θ_on, t_on = t_on);
        _δt = _δts[1];

        # run the budgets for all ∂x∂t
        θ_on ? soil_budget!(config, spac, _δt) : nothing;
        if spac._root_connection
            stomatal_conductance!(spac, _δt);
            t_on ? plant_energy!(config, spac, _δt) : nothing;
            p_on ? xylem_flow_profile!(spac, _δt) : nothing;
        end;

        _t_res -= _δt;

        # if _t_res > 0 rerun the budget functions (shortwave radiation not included) and etc., else break
        if _t_res > 0
            t_on ? longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL) : nothing;
            if spac._root_connection
                p_on ? xylem_pressure_profile!(spac; update = update) : nothing;
                leaf_photosynthesis!(spac, GCO₂Mode());
            end;
            θ_on ? soil_budget!(config, spac) : nothing;
            if spac._root_connection
                stomatal_conductance!(spac);
                t_on ? plant_energy!(config, spac) : nothing;
            end;
        else
            break;
        end;

        # if total count exceeds 100
        if (_count > 1000) && (_δt < 0.01) && (_t_res > 10)
            @info "Number of steppers exceeds 100, breaking..." spac.LATITUDE spac.LONGITUDE spac.CANOPY.lai _t_res _δts;
            break;
        end;
    end;

    return nothing
);

time_stepper!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}; update::Bool = false) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES, METEO, SOIL) = spac;

    # run the update function until the gpp is stable
    _count = 0;
    _gpp_last = -1;
    while true
        # compute the dxdt (not shortwave radiation simulation)
        longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL);
        if spac._root_connection
            xylem_pressure_profile!(spac; update = update);
            leaf_photosynthesis!(spac, GCO₂Mode());
        end;
        soil_budget!(config, spac);
        if spac._root_connection
            stomatal_conductance!(spac);
            plant_energy!(config, spac);
        end;

        # determine whether to break the while loop
        _gpp = GPP(spac);
        _count += 1;
        if abs(_gpp - _gpp_last) < 1e-6 || _count > 5000
            break;
        end;
        _gpp_last = _gpp;

        # use adjusted time to make sure no numerical issues
        _δts = adjusted_time(spac, config, FT(30); θ_on = false);
        _δt = _δts[1];

        # run the budgets for all ∂x∂t (except for soil)
        if spac._root_connection
            stomatal_conductance!(spac, _δt);
            plant_energy!(config, spac, _δt);
            xylem_flow_profile!(spac, _δt);
        end;
    end;

    # run canopy fluorescence
    canopy_fluorescence!(spac, config);

    return nothing
);

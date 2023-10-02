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
#     2023-Sep-11: move the optional t_on and θ_on to the config struct
#     2023-Sep-30: add time controller of junction water
#
#######################################################################################################################################################################################################
"""

        adjusted_time(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Return adjusted time that soil does not over saturate or drain, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
function adjusted_time(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    (; DEBUG, ENABLE_ENERGY_BUDGET, ENABLE_SOIL_WATER_BUDGET) = config;
    (; BRANCHES, JUNCTION, LEAVES, SOIL, TRUNK) = spac;

    _δt_1 = δt;

    # make sure each layer does not drain (allow for oversaturation), and θ change is less than 0.01
    _δt_2 = _δt_1;
    if ENABLE_SOIL_WATER_BUDGET
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
    if ENABLE_ENERGY_BUDGET
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
    if ENABLE_ENERGY_BUDGET
        _∂T∂t = TRUNK.energy.auxil.∂e∂t / (CP_L_MOL(FT) * sum(TRUNK.xylem.state.v_storage));
        _δt_4 = min(1 / abs(_∂T∂t), _δt_4);

        if DEBUG
            if any(isnan, (_δt_4, _∂T∂t))
                @info "Debugging" _δt_4 _∂T∂t;
                error("NaN in adjusted_time at trunk");
            end;
        end;
    end;

    # make sure the junction water does not change more than 10 mol per time step
    _δt_5 = _δt_4;
    if JUNCTION.auxil.∂w∂t != 0
        _δt_5 = min(10 / abs(JUNCTION.auxil.∂w∂t), _δt_5);

        if DEBUG
            if any(isnan, (_δt_5, JUNCTION.auxil.∂w∂t))
                @info "Debugging" _δt_5 JUNCTION.auxil.∂w∂t;
                error("NaN in adjusted_time at junction");
            end;
        end;
    end;

    # make sure branch stem temperatures do not change more than 1 K per time step
    _δt_6 = _δt_5;
    if ENABLE_ENERGY_BUDGET
        for _branch in BRANCHES
            _∂T∂t = _branch.energy.auxil.∂e∂t / (CP_L_MOL(FT) * sum(_branch.xylem.state.v_storage));
            _δt_6 = min(1 / abs(_∂T∂t), _δt_6);

            if DEBUG
                if any(isnan, (_δt_6, _∂T∂t))
                    @info "Debugging" _δt_6 _∂T∂t;
                    error("NaN in adjusted_time at branch");
                end;
            end;
        end;
    end;

    # make sure leaf temperatures do not change more than 1 K per time step
    _δt_7 = _δt_6;
    if ENABLE_ENERGY_BUDGET && spac.CANOPY.lai > 0
        for _clayer in LEAVES
            _∂T∂t = _clayer.NS.energy.auxil.∂e∂t / (_clayer.NS.xylem.state.cp * _clayer.NS.bio.state.lma * 10 + CP_L_MOL(FT) * _clayer.NS.capacitor.state.v_storage);
            _δt_7 = min(1 / abs(_∂T∂t), _δt_7);

            if DEBUG
                if any(isnan, (_δt_7, _∂T∂t))
                    @info "Debugging" _δt_7 _∂T∂t;
                    error("NaN in adjusted_time at leaf");
                end;
            end;
        end;
    end;

    # make sure leaf stomatal conductances do not change more than 0.06 mol m⁻² s⁻¹
    _δt_8 = _δt_7;
    if spac.CANOPY.lai > 0
        for _clayer in LEAVES
            for _∂g∂t in _clayer.∂g∂t_sunlit
                _δt_8 = min(FT(0.06) / abs(_∂g∂t), _δt_8);
            end;
            _δt_8 = min(FT(0.06) / abs(_clayer.∂g∂t_shaded), _δt_8);

            if DEBUG
                if any(isnan, (_δt_8, _clayer.∂g∂t_shaded, mean(_clayer.∂g∂t_sunlit)))
                    @info "Debugging" _δt_8 _clayer.∂g∂t_shaded mean(_clayer.∂g∂t_sunlit);
                    error("NaN in adjusted_time at leaf stomatal conductance");
                end;
            end;
        end;
    end;

    # make sure adjusted time is not nan
    @assert !isnan(_δt_8) "NaN in adjusted_time";

    return (_δt_8, _δt_7, _δt_6, _δt_5, _δt_4, _δt_3, _δt_2, _δt_1)
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
#     2023-Sep-11: move the optional p_on, t_on and θ_on to the config struct
#     2023-Sep-14: remove some if else control from root disconnection
#
#######################################################################################################################################################################################################
"""

    time_stepper!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::Number) where {FT}
    time_stepper!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Move forward in time for SPAC with time stepper controller, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step (if not given, solve for steady state solution)

"""
function time_stepper! end

time_stepper!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::Number) where {FT} = (
    (; CANOPY, LEAVES, METEO, SOIL) = spac;

    # run the update function until time elapses
    _count = 0;
    _t_res = FT(δt);
    while true
        _count += 1;
        _δts = adjusted_time(config, spac, _t_res);
        _δt = _δts[1];

        # run the budgets for all ∂x∂t
        soil_budget!(config, spac, _δt);
        stomatal_conductance!(spac, _δt);
        plant_energy!(config, spac, _δt);
        if spac._root_connection
            plant_water_budget!(spac, _δt);
        end;

        _t_res -= _δt;

        # if _t_res > 0 rerun the budget functions (shortwave radiation not included) and etc., else break
        if _t_res > 0
            update_substep_auxils!(spac);
            longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL);
            if spac._root_connection
                plant_flow_profile!(config, spac);
                plant_pressure_profile!(config, spac);
            end;
            leaf_photosynthesis!(spac, GCO₂Mode());
            soil_budget!(config, spac);
            stomatal_conductance!(spac);
            spac_energy_flow!(spac);
        else
            break;
        end;

        # if total count exceeds 100
        if (_count > 1000) && (_δt < 0.01) && (_t_res > 10)
            @info "Number of steppers exceeds 1000, breaking..." spac.LATITUDE spac.LONGITUDE spac.CANOPY.lai _t_res _δts;
            break;
        end;
    end;

    return nothing
);

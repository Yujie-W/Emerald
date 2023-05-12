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
#
#######################################################################################################################################################################################################
"""

    function adjusted_time(spac::MultiLayerSPAC{FT}, δt::FT; t_on::Bool = true, θ_on::Bool = true) where {FT<:AbstractFloat}

Return adjusted time that soil does not over saturate or drain, given
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step
- `t_on` If true, plant energy budget is on (set false to run sensitivity analysis or prescribing mode)
- `θ_on` If true, soil water budget is on (set false to run sensitivity analysis or prescribing mode)

"""
function adjusted_time(spac::MultiLayerSPAC{FT}, δt::FT; t_on::Bool = true, θ_on::Bool = true) where {FT<:AbstractFloat}
    (; BRANCHES, DIM_LAYER, LEAVES, SOIL, TRUNK) = spac;

    _δt_1 = δt;

    # make sure each layer does not drain (allow for oversaturation), and θ change is less than 0.01
    _δt_2 = _δt_1;
    if θ_on
        for _i in 1:SOIL.DIM_SOIL
            _δt_2 = min(FT(0.01) / abs(SOIL.LAYERS[_i].∂θ∂t), _δt_2);
            if SOIL.LAYERS[_i].∂θ∂t < 0
                _δt_dra = (SOIL.LAYERS[_i].VC.Θ_RES - SOIL.LAYERS[_i].θ) / SOIL.LAYERS[_i].∂θ∂t;
                _δt_2 = min(_δt_dra, _δt_2);
            end;
        end;
    end;

    # make sure soil temperatures do not change more than 1 K per time step
    _δt_3 = _δt_2;
    if t_on
        for _i in 1:SOIL.DIM_SOIL
            _∂T∂t = SOIL.LAYERS[_i].∂e∂t / (SOIL.LAYERS[_i].CP * SOIL.LAYERS[_i].ρ + SOIL.LAYERS[_i].θ * ρ_H₂O(FT) * CP_L(FT));
            _δt_3 = min(1 / abs(_∂T∂t), _δt_3);
        end;
    end;

    # make sure trunk temperatures do not change more than 1 K per time step
    _δt_4 = _δt_3;
    if t_on
        _∂T∂t = TRUNK.∂e∂t / (CP_L_MOL(FT) * sum(TRUNK.HS.v_storage));
        _δt_4 = min(1 / abs(_∂T∂t), _δt_4);
    end;

    # make sure branch stem temperatures do not change more than 1 K per time step
    _δt_5 = _δt_4;
    if t_on
        for _i in 1:DIM_LAYER
            _∂T∂t = BRANCHES[_i].∂e∂t / (CP_L_MOL(FT) * sum(BRANCHES[_i].HS.v_storage));
            _δt_5 = min(1 / abs(_∂T∂t), _δt_5);
        end;
    end;

    # make sure leaf temperatures do not change more than 1 K per time step
    _δt_6 = _δt_5;
    if t_on
        for _i in 1:DIM_LAYER
            _∂T∂t = LEAVES[_i].∂e∂t / (LEAVES[_i].CP * LEAVES[_i].BIO.lma * 10 + CP_L_MOL(FT) * LEAVES[_i].HS.v_storage);
            _δt_6 = min(1 / abs(_∂T∂t), _δt_6);
        end;
    end;

    # make sure leaf stomatal conductances do not change more than 0.01 mol m⁻² s⁻¹
    _δt_7 = _δt_6;
    for _i in 1:DIM_LAYER
        for _∂g∂t in LEAVES[_i].∂g∂t_sunlit
            _δt_7 = min(FT(0.06) / abs(_∂g∂t), _δt_7);
        end;
        _δt_7 = min(FT(0.06) / abs(LEAVES[_i].∂g∂t_shaded), _δt_7);
    end;

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
#
#######################################################################################################################################################################################################
"""

    time_stepper!(spac::MultiLayerSPAC{FT}, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat}
    time_stepper!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}; update::Bool = false) where {FT<:AbstractFloat}

Move forward in time for SPAC with time stepper controller, given
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step (if not given, solve for steady state solution)
- `p_on` If true, plant hydraulic flow and pressure profiles will be updated
- `t_on` If true, plant energy budget is on (set false to run sensitivity analysis or prescribing mode)
- `update` If true, update leaf xylem legacy effect
- `θ_on` If true, soil water budget is on (set false to run sensitivity analysis or prescribing mode)
- `config` Configuration for `MultiLayerSPAC`

"""
function time_stepper! end

time_stepper!(spac::MultiLayerSPAC{FT}, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES, METEO, SOIL) = spac;

    # run the update function until time elapses
    _count = 0;
    _t_res = FT(δt);
    while true
        _count += 1;
        _δts = adjusted_time(spac, _t_res; θ_on = θ_on, t_on = t_on);
        _δt = _δts[1];

        # run the budgets for all ∂x∂t
        if θ_on soil_budget!(spac, _δt); end;
        stomatal_conductance!(spac, _δt);
        if t_on plant_energy!(spac, _δt); end;
        if p_on xylem_flow_profile!(spac, _δt); end;

        _t_res -= _δt;

        # if _t_res > 0 rerun the budget functions (shortwave radiation not included) and etc., else break
        if _t_res > 0
            if t_on longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL); end;
            if p_on xylem_pressure_profile!(spac; update = update); end;
            leaf_photosynthesis!(spac, GCO₂Mode());
            if θ_on soil_budget!(spac); end;
            stomatal_conductance!(spac);
            if t_on plant_energy!(spac); end;
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
        xylem_pressure_profile!(spac; update = update);
        leaf_photosynthesis!(spac, GCO₂Mode());
        soil_budget!(spac);
        stomatal_conductance!(spac);
        plant_energy!(spac);

        # determine whether to break the while loop
        _gpp = GPP(spac);
        _count += 1;
        if abs(_gpp - _gpp_last) < 1e-6 || _count > 5000
            break;
        end;
        _gpp_last = _gpp;

        # use adjusted time to make sure no numerical issues
        _δts = adjusted_time(spac, FT(30); θ_on = false);
        _δt = _δts[1];

        # run the budgets for all ∂x∂t (except for soil)
        stomatal_conductance!(spac, _δt);
        plant_energy!(spac, _δt);
        xylem_flow_profile!(spac, _δt);
    end;

    # run canopy fluorescence
    canopy_fluorescence!(spac, config);

    return nothing
);
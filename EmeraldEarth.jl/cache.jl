#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: initialize CACHE_STATE at the same time
#     2023-Mar-13: initialize CACHE_CONFIG at the same time
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
# Bug fixes
#     2023-Aug-26: make sure sza < 89 when total radiation is higher than 10 W m⁻²
#
#######################################################################################################################################################################################################
"""

    initialize_cache!(FT)

Initialize the global parameter `CACHE_SPAC`, given
- `FT` Floating type

"""
function initialize_cache!(FT)
    global CACHE_CONFIG, CACHE_SPAC, CACHE_STATE;

    # create a SPAC to work on
    _z_canopy = FT(10);
    CACHE_CONFIG = SPACConfiguration{FT}();
    CACHE_SPAC = MultiLayerSPAC(
                CACHE_CONFIG;
                air_bounds = collect(0:21) * _z_canopy / 20,
                latitude = 0,
                longitude = 0,
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                zs = [-2, _z_canopy/2, _z_canopy]);

    # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
    # for _organ in [CACHE_SPAC.LEAVES; CACHE_SPAC.BRANCHES; CACHE_SPAC.TRUNK; CACHE_SPAC.ROOTS]
    #     _organ.HS.VC.B = 3;
    #     _organ.HS.VC.C = 1;
    # end;

    # update leaf mass per area and stomtal model
    @inline linear_p_soil(x) = min(1, max(eps(FT), 1 + x / 5));
    _bt = BetaFunction{FT}(FUNC = linear_p_soil, PARAM_X = BetaParameterPsoil(), PARAM_Y = BetaParameterG1());
    for _leaves in CACHE_SPAC.LEAVES
        _leaves.SM = MedlynSM{FT}(G0 = 0.005, β = _bt);
    end;

    # initialize the spac with non-saturated soil
    update!(CACHE_SPAC, CACHE_CONFIG; swcs = Tuple(max(_slayer.VC.Θ_SAT - 0.02, (_slayer.VC.Θ_SAT + _slayer.VC.Θ_RES) / 2) for _slayer in CACHE_SPAC.SOIL.LAYERS));
    initialize!(CACHE_SPAC, CACHE_CONFIG);

    # create a state struct based on the spac
    CACHE_STATE = MultiLayerSPACState{FT}(CACHE_SPAC);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: add step to synchronize state variables into CACHE_SPAC
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Apr-13: re-wire RAD_SW_REF to CACHE_CONFIG
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-15: make sure prescribed soil parameters are not NaN and rad is >= 0
#
#######################################################################################################################################################################################################
"""

    synchronize_cache!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MultiLayerSPACState{FT}}) where {FT<:AbstractFloat}

Synchronize SPAC parameters from,
- `gm_params` Dict for GriddingMachine parameters
- `wd_params` Dict for weather drivers
- `state` `MultiLayerSPACState` for all state variables, or nothing

"""
function synchronize_cache!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MultiLayerSPACState})
    FT = gm_params["FT"];
    _z_canopy = max(FT(0.1), gm_params["CANOPY_HEIGHT"]);

    #
    # TODO: update canopy height for plant hydraulic system based on _z_canopy
    #

    # update the values in the CACHE_SPAC, use .= for arrays
    global CACHE_SPAC;
    CACHE_SPAC.LATITUDE = gm_params["LATITUDE"];
    CACHE_SPAC.LONGITUDE = gm_params["LONGITUDE"];
    CACHE_SPAC.ELEVATION = gm_params["ELEVATION"];
    CACHE_SPAC.Z .= [-2, _z_canopy/2, _z_canopy];
    CACHE_SPAC.Z_AIR .= collect(0:21) * _z_canopy / 20;
    CACHE_SPAC.SOIL.COLOR = gm_params["SOIL_COLOR"];

    # update soil type information per layer
    for _i in eachindex(CACHE_SPAC.SOIL.LAYERS)
        # TODO: add a line to parameterize K_MAX
        # TODO: fix these later with better data source
        if !isnan(gm_params["SOIL_α"][_i]) && !isnan(gm_params["SOIL_N"][_i]) && !isnan(gm_params["SOIL_ΘR"][_i]) && !isnan(gm_params["SOIL_ΘS"][_i])
            CACHE_SPAC.SOIL.LAYERS[_i].VC.α = gm_params["SOIL_α"][_i];
            CACHE_SPAC.SOIL.LAYERS[_i].VC.N = gm_params["SOIL_N"][_i];
            CACHE_SPAC.SOIL.LAYERS[_i].VC.M = 1 - 1 / CACHE_SPAC.SOIL.LAYERS[_i].VC.N;
            CACHE_SPAC.SOIL.LAYERS[_i].VC.Θ_RES = gm_params["SOIL_ΘR"][_i];
            CACHE_SPAC.SOIL.LAYERS[_i].VC.Θ_SAT = gm_params["SOIL_ΘS"][_i];
        end;
    end;

    # update leaf mass per area and stomtal model
    for _leaves in CACHE_SPAC.LEAVES
        _leaves.BIO.lma = gm_params["LMA"];
        _leaves.SM.G1 = gm_params["MEDLYN_G1"];
    end;

    # update environmental conditions
    for _alayer in CACHE_SPAC.AIR
        _alayer.P_AIR = wd_params["P_ATM"];
        update!(_alayer; t = wd_params["T_AIR"], vpd = wd_params["VPD"], wind = wd_params["WIND"]);
    end;

    # sync the environmental conditions per layer for CO₂ concentration
    if !isnothing(gm_params["CO2"])
        for _alayer in CACHE_SPAC.AIR
            update!(_alayer; f_CO₂ = gm_params["CO2"]);
        end;
    end;

    # update shortwave and longwave radiation
    _in_dir = CACHE_CONFIG.RAD_SW_REF.e_direct'  * CACHE_CONFIG.WLSET.ΔΛ / 1000;
    _in_dif = CACHE_CONFIG.RAD_SW_REF.e_diffuse' * CACHE_CONFIG.WLSET.ΔΛ / 1000;
    CACHE_SPAC.METEO.rad_sw.e_direct  .= CACHE_CONFIG.RAD_SW_REF.e_direct  .* max(0,wd_params["RAD_DIR"]) ./ _in_dir;
    CACHE_SPAC.METEO.rad_sw.e_diffuse .= CACHE_CONFIG.RAD_SW_REF.e_diffuse .* max(0,wd_params["RAD_DIF"]) ./ _in_dif;
    CACHE_SPAC.METEO.rad_lw = wd_params["RAD_LW"];

    # update solar zenith angle based on the time
    _sza = solar_zenith_angle(CACHE_SPAC.LATITUDE, FT(wd_params["FDOY"]));
    CACHE_SPAC.ANGLES.sza = (wd_params["RAD_DIR"] + wd_params["RAD_DIF"] > 10) ? min(_sza, 88.999) : _sza;

    # prescribe soil water content
    if "SWC" in keys(wd_params)
        update!(CACHE_SPAC, CACHE_CONFIG; swcs = wd_params["SWC"], t_soils = wd_params["T_SOIL"]);
        # TODO: remove this when soil diffusion problem is fixed
        update!(CACHE_SPAC, CACHE_CONFIG; swcs = Tuple(min(_slayer.VC.Θ_SAT - 0.001, _slayer.θ) for _slayer in CACHE_SPAC.SOIL.LAYERS));
    end;

    # synchronize the state if state is not nothing, otherwise set all values to NaN (do thing before prescribing T_SKIN)
    if !isnothing(state)
        spac_state!(state, CACHE_SPAC);
    else
        CACHE_SPAC.MEMORY.tem .= NaN;
    end;

    # prescribe leaf temperature from skin temperature
    # TODO: add CACHE_SPAC.MEMORY.tem as prognostic variable
    if "T_SKIN" in keys(wd_params)
        push!(CACHE_SPAC.MEMORY.tem, wd_params["T_SKIN"]);
        if length(CACHE_SPAC.MEMORY.tem) > 240 deleteat!(CACHE_SPAC.MEMORY.tem,1) end;
        update!(CACHE_SPAC, CACHE_CONFIG; t_leaf = wd_params["T_SKIN"], t_clm = nanmean(CACHE_SPAC.MEMORY.tem));
    end;

    # synchronize LAI, CHL, and CI
    _iday = Int(floor(wd_params["INDEX"] / 24)) + 1;
    _chl = griddingmachine_data(gm_params["CHLOROPHYLL"], gm_params["YEAR"], _iday);
    _cli = griddingmachine_data(gm_params["CLUMPING"], gm_params["YEAR"], _iday);
    _lai = griddingmachine_data(gm_params["LAI"], gm_params["YEAR"], _iday);
    _vcm = griddingmachine_data(gm_params["VCMAX25"], gm_params["YEAR"], _iday);

    # update clumping index, LAI, Vcmax, and Chl
    update!(CACHE_SPAC, CACHE_CONFIG; cab = _chl, car = _chl / 7, ci = _cli, lai =_lai, vcmax = _vcm, vcmax_expo = 0.3);

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to determine the day bounds of input GriddingMachine drivers
# TODO:
#     Generalize this function in EmeraldMath or EmeraldUtility
#
#######################################################################################################################################################################################################
"""

    griddingmachine_data(data::Vector, year::Int, d::Int)

Return the index of data, given
- `data` Time series of input data
- `year` Year
- `d` Day number

"""
function griddingmachine_data(data::Vector, year::Int, d::Int)
    _bounds = [0,367];
    _n = length(data);
    if _n == 1
        _bounds = [0,367]
    elseif _n == 12
        _bounds = isleapyear(year) ? MDAYS_LEAP : MDAYS;
    elseif _n == 46
        _bounds = [collect(0:8:361); 367]
    elseif _n == 52
        _bounds = [collect(0:7:361); 367]
    else
        error("This temporal resolution is not supported: $(_n)!");
    end

    _ind = findfirst(d .<= _bounds) - 1

    return data[_ind]
end

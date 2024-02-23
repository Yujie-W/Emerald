


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: add step to synchronize state variables into CACHE_SPAC
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-15: make sure prescribed soil parameters are not NaN and rad is >= 0
#
#######################################################################################################################################################################################################
"""

    synchronize_cache!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MultiLayerSPACState{FT}})

Synchronize SPAC parameters from,
- `gm_params` Dict for GriddingMachine parameters
- `wd_params` Dict for weather drivers
- `state` `MultiLayerSPACState` for all state variables, or nothing

"""
function synchronize_cache!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing})
    FT = gm_params["FT"];
    _z_canopy = max(FT(0.1), gm_params["CANOPY_HEIGHT"]);

    #
    # TODO: update canopy height for plant hydraulic system based on _z_canopy
    #

    # update the values in the CACHE_SPAC, use .= for arrays
    global CACHE_SPAC;
    CACHE_SPAC.info.lat = gm_params["LATITUDE"];
    CACHE_SPAC.info.lon = gm_params["LONGITUDE"];
    CACHE_SPAC.info.elev = gm_params["ELEVATION"];
    CACHE_SPAC.Z .= [-2, _z_canopy/2, _z_canopy];
    CACHE_SPAC.Z_AIR .= collect(0:21) * _z_canopy / 20;
    CACHE_SPAC.soil_bulk.state.color = gm_params["SOIL_COLOR"];

    # update soil type information per layer
    for i in eachindex(CACHE_SPAC.soils)
        # TODO: add a line to parameterize K_MAX
        # TODO: fix these later with better data source
        if !isnan(gm_params["SOIL_α"][i]) && !isnan(gm_params["SOIL_N"][i]) && !isnan(gm_params["SOIL_ΘR"][i]) && !isnan(gm_params["SOIL_ΘS"][i])
            CACHE_SPAC.soils[i].state.vc.α = gm_params["SOIL_α"][i];
            CACHE_SPAC.soils[i].state.vc.N = gm_params["SOIL_N"][i];
            CACHE_SPAC.soils[i].state.vc.M = 1 - 1 / CACHE_SPAC.soils[i].state.vc.N;
            CACHE_SPAC.soils[i].state.vc.Θ_RES = gm_params["SOIL_ΘR"][i];
            CACHE_SPAC.soils[i].state.vc.Θ_SAT = gm_params["SOIL_ΘS"][i];
        end;
    end;

    # update leaf mass per area and stomtal model
    for leaf in CACHE_SPAC.plant.leaves
        leaf.BIO.state.lma = gm_params["LMA"];
        leaf.flux.state.stomatal_model.G1 = gm_params["MEDLYN_G1"];
    end;

    # update environmental conditions
    for air in CACHE_SPAC.airs
        air.state.p_air = wd_params["P_ATM"];
        prescribe_air!(air; t = wd_params["T_AIR"], vpd = wd_params["VPD"], wind = wd_params["WIND"]);
    end;

    # sync the environmental conditions per layer for CO₂ concentration
    if !isnothing(gm_params["CO2"])
        for air in CACHE_SPAC.airs
            prescribe_air!(air; f_CO₂ = gm_params["CO2"]);
        end;
    end;

    # update shortwave and longwave radiation
    _in_dir = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1)'  * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    _in_dif = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2)' * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    CACHE_SPAC.meteo.rad_sw.e_dir .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1) .* max(0,wd_params["RAD_DIR"]) ./ _in_dir;
    CACHE_SPAC.meteo.rad_sw.e_dif .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2) .* max(0,wd_params["RAD_DIF"]) ./ _in_dif;
    CACHE_SPAC.meteo.rad_lw = wd_params["RAD_LW"];

    # update solar zenith angle based on the time




    _sza = solar_zenith_angle(CACHE_SPAC.info.lat, FT(wd_params["FDOY"]));




    CACHE_SPAC.canopy.sun_geometry.state.sza = (wd_params["RAD_DIR"] + wd_params["RAD_DIF"] > 10) ? min(_sza, 88.999) : _sza;



    # prescribe soil water content
    if "SWC" in keys(wd_params)
        prescribe_soil!(CACHE_SPAC; swcs = wd_params["SWC"], t_soils = wd_params["T_SOIL"]);
        # TODO: remove this when soil diffusion problem is fixed
        prescribe_soil!(CACHE_SPAC; swcs = Tuple(min(soil.state.vc.Θ_SAT - 0.001, soil.state.θ) for soil in CACHE_SPAC.soils));
    end;

    # synchronize the state if state is not nothing, otherwise set all values to NaN (do thing before prescribing T_SKIN)
    if !isnothing(state)
        #spac_state!(state, CACHE_SPAC);
    else
        CACHE_SPAC.plant.memory.t_history .= NaN;
    end;

    # prescribe leaf temperature from skin temperature
    if "T_SKIN" in keys(wd_params)
        push!(CACHE_SPAC.plant.memory.t_history, wd_params["T_SKIN"]);
        if length(CACHE_SPAC.plant.memory.t_history) > 240 deleteat!(CACHE_SPAC.plant.memory.t_history,1) end;
        prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; t_leaf = wd_params["T_SKIN"], t_clm = nanmean(CACHE_SPAC.plant.memory.t_history));
    end;

    # synchronize LAI, CHL, and CI
    _iday = Int(floor(wd_params["INDEX"] / 24)) + 1;
    _chl = griddingmachine_data(gm_params["CHLOROPHYLL"], gm_params["YEAR"], _iday);
    _cli = griddingmachine_data(gm_params["CLUMPING"], gm_params["YEAR"], _iday);
    _lai = griddingmachine_data(gm_params["LAI"], gm_params["YEAR"], _iday);
    _vcm = griddingmachine_data(gm_params["VCMAX25"], gm_params["YEAR"], _iday);

    # update clumping index, LAI, Vcmax, and Chl
    prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; cab = _chl, car = _chl / 7, ci = _cli, lai =_lai, vcmax = _vcm, vcmax_expo = 0.3);

    return nothing
end;


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
    end;

    _ind = findfirst(d .<= _bounds) - 1

    return data[_ind]
end;

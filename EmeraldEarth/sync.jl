


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

    synchronize_cache!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::Union{Nothing,MultiLayerSPACState{FT}})

Synchronize SPAC parameters from,
- `gm_dict` Dict for GriddingMachine parameters
- `wd_dict` Dict for weather drivers
- `state` `MultiLayerSPACState` for all state variables, or nothing

"""
function synchronize_cache!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::Union{Nothing})
    FT = gm_dict["FT"];
    _z_canopy = max(FT(0.1), gm_dict["CANOPY_HEIGHT"]);

    #
    # TODO: update canopy height for plant hydraulic system based on _z_canopy
    #

    # update the values in the CACHE_SPAC, use .= for arrays
    global CACHE_SPAC;
    CACHE_SPAC.info.lat = gm_dict["LATITUDE"];
    CACHE_SPAC.info.lon = gm_dict["LONGITUDE"];
    CACHE_SPAC.info.elev = gm_dict["ELEVATION"];
    CACHE_SPAC.Z .= [-2, _z_canopy/2, _z_canopy];
    CACHE_SPAC.Z_AIR .= collect(0:21) * _z_canopy / 20;
    CACHE_SPAC.soil_bulk.state.color = gm_dict["SOIL_COLOR"];

    # update soil type information per layer
    for i in eachindex(CACHE_SPAC.soils)
        # TODO: add a line to parameterize K_MAX
        # TODO: fix these later with better data source
        if !isnan(gm_dict["SOIL_α"][i]) && !isnan(gm_dict["SOIL_N"][i]) && !isnan(gm_dict["SOIL_ΘR"][i]) && !isnan(gm_dict["SOIL_ΘS"][i])
            CACHE_SPAC.soils[i].state.vc.α = gm_dict["SOIL_α"][i];
            CACHE_SPAC.soils[i].state.vc.N = gm_dict["SOIL_N"][i];
            CACHE_SPAC.soils[i].state.vc.M = 1 - 1 / CACHE_SPAC.soils[i].state.vc.N;
            CACHE_SPAC.soils[i].state.vc.Θ_RES = gm_dict["SOIL_ΘR"][i];
            CACHE_SPAC.soils[i].state.vc.Θ_SAT = gm_dict["SOIL_ΘS"][i];
        end;
    end;

    # update leaf mass per area and stomtal model
    for leaf in CACHE_SPAC.plant.leaves
        leaf.BIO.state.lma = gm_dict["LMA"];
        leaf.flux.state.stomatal_model.G1 = gm_dict["G1_MEDLYN_C3"];
    end;

    # update environmental conditions
    for air in CACHE_SPAC.airs
        air.state.p_air = wd_dict["P_ATM"];
        prescribe_air!(air; t = wd_dict["T_AIR"], vpd = wd_dict["VPD"], wind = wd_dict["WIND"]);
    end;

    # sync the environmental conditions per layer for CO₂ concentration
    if !isnothing(gm_dict["CO2"])
        for air in CACHE_SPAC.airs
            prescribe_air!(air; f_CO₂ = gm_dict["CO2"]);
        end;
    end;

    # update shortwave and longwave radiation
    _in_dir = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1)'  * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    _in_dif = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2)' * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    CACHE_SPAC.meteo.rad_sw.e_dir .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1) .* max(0,wd_dict["RAD_DIR"]) ./ _in_dir;
    CACHE_SPAC.meteo.rad_sw.e_dif .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2) .* max(0,wd_dict["RAD_DIF"]) ./ _in_dif;
    CACHE_SPAC.meteo.rad_lw = wd_dict["RAD_LW"];

    # update solar zenith angle based on the time




    _sza = solar_zenith_angle(CACHE_SPAC.info.lat, FT(wd_dict["FDOY"]));




    CACHE_SPAC.canopy.sun_geometry.state.sza = (wd_dict["RAD_DIR"] + wd_dict["RAD_DIF"] > 10) ? min(_sza, 88.999) : _sza;



    # prescribe soil water content
    if "SWC" in keys(wd_dict)
        prescribe_soil!(CACHE_SPAC; swcs = wd_dict["SWC"], t_soils = wd_dict["T_SOIL"]);
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
    if "T_SKIN" in keys(wd_dict)
        push!(CACHE_SPAC.plant.memory.t_history, wd_dict["T_SKIN"]);
        if length(CACHE_SPAC.plant.memory.t_history) > 240 deleteat!(CACHE_SPAC.plant.memory.t_history,1) end;
        prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; t_leaf = wd_dict["T_SKIN"], t_clm = nanmean(CACHE_SPAC.plant.memory.t_history));
    end;

    # synchronize LAI, CHL, and CI
    _iday = Int(floor(wd_dict["INDEX"] / 24)) + 1;
    _chl = query_griddingmachine_data(gm_dict["CHLOROPHYLL"], gm_dict["YEAR"], _iday);
    _cli = query_griddingmachine_data(gm_dict["CLUMPING"], gm_dict["YEAR"], _iday);
    _lai = query_griddingmachine_data(gm_dict["LAI"], gm_dict["YEAR"], _iday);
    _vcm = query_griddingmachine_data(gm_dict["VCMAX25"], gm_dict["YEAR"], _iday);

    # update clumping index, LAI, Vcmax, and Chl
    prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; ci = _cli, lai =_lai, vcmax = _vcm, vcmax_expo = 0.3);
    prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; cab = _chl, car = _chl / 7);

    return nothing
end;

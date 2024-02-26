#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to prescribe parameters from weather drivers
#     2023-Mar-27: prescribe T only if t_on is true, prescribe SWC only is θ_on is true
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: add a controller over rad to make sure it is >= 0
#     2023-Aug-25: make soil and leaf temperature and soil moisture optional
#     2024-Feb-22: remove state and auxil from spac struct
# Bug fixes
#     2023-Aug-26: computed sza in the middle of a time pierod may be > 0 when cumulated radiation is > 0, set it to 88.999
#
#######################################################################################################################################################################################################
"""

    prescribe!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, dfr::DataFrameRow; t_on::Bool = true, θ_on::Bool = true) where {FT}

Prescribe traits and environmental conditions, given
- `config` `SPACConfiguration` type SPAC configuration
- `spac` `BulkSPAC` type SPAC
- `dfr` `DataFrameRow` type weather driver

"""
function prescribe!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, dfr::DataFrameRow; initialize_state::Bool = false) where {FT}
    # read the data out of dataframe row to reduce memory allocation
    _df_atm::FT = dfr.P_ATM;
    _df_chl::FT = dfr.CHLOROPHYLL;
    _df_cli::FT = dfr.CI;
    _df_co2::FT = dfr.CO2;
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;
    _df_doy::FT = dfr.FDOY;
    _df_lai::FT = dfr.LAI;
    _df_lwr::FT = dfr.RAD_LW;
    _df_pcp::FT = dfr.PRECIP;
    _df_tar::FT = dfr.T_AIR;
    _df_vcm::FT = dfr.VCMAX25;
    _df_vpd::FT = dfr.VPD;
    _df_wnd::FT = dfr.WIND;

    # adjust optimum t based on 10 day moving average skin temperature
    prescribe_traits!(config, spac; t_clm = nanmean(spac.plant.memory.t_history));

    # prescribe soil water contents and leaf temperature (for version B1 only)
    if initialize_state
        _df_tlf::FT = dfr.T_LEAF;
        _df_ts1::FT = dfr.T_SOIL_1;
        _df_ts2::FT = dfr.T_SOIL_2;
        _df_ts3::FT = dfr.T_SOIL_3;
        _df_ts4::FT = dfr.T_SOIL_4;
        prescribe_traits!(config, spac; t_leaf = max(_df_tar, _df_tlf));

        _df_sw1::FT = dfr.SWC_1;
        _df_sw2::FT = dfr.SWC_2;
        _df_sw3::FT = dfr.SWC_3;
        _df_sw4::FT = dfr.SWC_4;
        prescribe_soil!(spac; swcs = (_df_sw1, _df_sw2, _df_sw3, _df_sw4), t_soils = (_df_ts1, _df_ts2, _df_ts3, _df_ts4));
    end;

    # prescribe the precipitation related parameters
    spac.meteo.rain = _df_pcp * ρ_H₂O(FT) / M_H₂O(FT) / 3600;
    spac.meteo.t_precip = _df_tar;

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    _trigger_lai::Bool = !isnan(_df_lai) && (_df_lai != spac.plant.memory.lai);
    _trigger_vcm::Bool = !isnan(_df_vcm) && (_df_vcm != spac.plant.memory.vcmax25);
    _trigger_chl::Bool = !isnan(_df_chl) && (_df_chl != spac.plant.memory.chl);
    if _trigger_lai
        prescribe_traits!(config, spac; lai = _df_lai, vcmax_expo = 0.3);
        spac.plant.memory.lai = _df_lai;
    end;

    if _trigger_vcm
        prescribe_traits!(config, spac; vcmax = _df_vcm, vcmax_expo = 0.3);
        spac.plant.memory.vcmax25 = _df_vcm;
    end;

    if _trigger_chl
        prescribe_traits!(config, spac; cab = _df_chl, car = _df_chl / 7);
        spac.plant.memory.chl = _df_chl;
    end;

    # update clumping index
    prescribe_traits!(config, spac; ci = _df_cli);

    # update environmental conditions
    for air in spac.airs
        air.state.p_air = _df_atm;
        prescribe_air!(air; f_CO₂ = _df_co2, t = _df_tar, vpd = _df_vpd, wind = _df_wnd);
    end;

    # update downward shortwave and longwave radiation
    _in_dir = view(config.SPECTRA.SOLAR_RAD,:,1)' * config.SPECTRA.ΔΛ / 1000;
    _in_dif = view(config.SPECTRA.SOLAR_RAD,:,2)' * config.SPECTRA.ΔΛ / 1000;
    spac.meteo.rad_sw.e_dir .= view(config.SPECTRA.SOLAR_RAD,:,1) .* max(0,_df_dir) ./ _in_dir;
    spac.meteo.rad_sw.e_dif .= view(config.SPECTRA.SOLAR_RAD,:,2) .* max(0,_df_dif) ./ _in_dif;
    spac.meteo.rad_lw = _df_lwr;

    # update solar zenith angle based on the time
    _sza = solar_zenith_angle(spac.info.lat, FT(_df_doy));
    spac.canopy.sun_geometry.state.sza = (_df_dir + _df_dif > 10) ? min(_sza, 88.999) : _sza;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: move function from ClimaLand-0.2
#     2023-Mar-25: set reflectance based value to NaN at night
#     2023-Mar-27: add p_on, t_on, and θ_on options as in spac! function
#     2023-Mar-28: if option saving is false, return the simulated result dataframe
#     2023-Mar-28: add option selection to run part of the whole year simulations
#     2023-Mar-28: save swcs and temperatures based on t_on and θ_on
#     2023-Mar-29: add option to load initialial state from weather driver
#     2023-Aug-25: add method to run spac simulations using externally prepared variables
#     2023-Aug-26: add debug information
#     2023-Aug-27: show ind at debug mode, otherwise show progress bar
#     2023-Sep-07: initialize integrators when starting a new simulation in a long time step
#     2023-Sep-09: save the quantum yields when saving the simulation results
#     2023-Sep-11: save the integrated SIF when saving the simulation results
#
#######################################################################################################################################################################################################
"""

    simulation!(wd_tag::String, gmdict::Dict{String,Any}; appending::Bool = false, initialize_state::Union{Nothing,Bool} = true, saving::Union{Nothing,String} = nothing, selection = :)
    simulation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, wdf::DataFrame; initialize_state::Union{Nothing,Bool} = true, saving::Union{Nothing,String} = nothing, selection = :) where {FT}

Run simulation on site level, given
- `wd_tag` Weather drive tag such as `wd1`
- `gmdict` GriddingMachine dict for site information
- `appending` If true, append new variables to weather driver when querying the file (set it to true when encountering any errors)
- `initialize_state` Initial state of spac: if is a bool, load the first data from the weather driver
- `saving` If is not nothing, save the simulations as a Netcdf file in the working directory; if is nothing, return the simulated result dataframe
- `selection` Run selection of data, default is : (namely 1:end;)

The second method can be used to run externally prepared config, spac, and weather driver, given
- `config` SPAC configuration
- `spac` SPAC
- `wdf` Weather driver dataframe

"""
function simulation! end;

simulation!(wd_tag::String, gmdict::Dict{String,Any}; appending::Bool = false, initialize_state::Union{Nothing,Bool} = true, saving::Union{Nothing,String} = nothing, selection = :) = (
    config = spac_config(gmdict);
    spac = spac_struct(gmdict, config);
    wdf = weather_driver(wd_tag, gmdict; appending = appending);

    simulation!(config, spac, wdf; initialize_state = initialize_state, saving = saving, selection = selection);

    return isnothing(saving) ? wdf : nothing
);

simulation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, wdf::DataFrame; initialize_state::Union{Nothing,Bool} = true, saving::Union{Nothing,String} = nothing, selection = :) where {FT} = (
    (; MESSAGE_LEVEL) = config;

    _wdfr = eachrow(wdf);

    # initialize spac based on initialize_state
    prescribe!(config, spac, _wdfr[1]; initialize_state = initialize_state);

    # iterate through the time steps
    if MESSAGE_LEVEL == 0
        for _dfr in _wdfr[selection]
            simulation!(config, spac, _dfr);
        end;
    elseif MESSAGE_LEVEL == 1
        @showprogress for _dfr in _wdfr[selection]
            simulation!(config, spac, _dfr);
        end;
    elseif MESSAGE_LEVEL == 2
        for _dfr in _wdfr[selection]
            @show _dfr.ind;
            @time simulation!(config, spac, _dfr);
        end;
    else
        error("MESSAGE_LEVEL should be 0, 1, or 2");
    end;

    # save simulation results to hard drive
    if !isnothing(saving)
        save_nc!(saving, wdf[selection,[_n != "ind" for _n in names(wdf)]]);
    end;

    return nothing
);

simulation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, dfr::DataFrameRow; δt::Number = 3600) where {FT} = (
    # read the data out of dataframe row to reduce memory allocation
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;

    # prescribe parameters
    prescribe!(config, spac, dfr);

    # run the model
    soil_plant_air_continuum!(config, spac, δt);
    mean_tleaf = nanmean([l.energy.s_aux.t for l in spac.plant.leaves]);
    push!(spac.plant.memory.t_history, mean_tleaf);
    if length(spac.plant.memory.t_history) > 240 deleteat!(spac.plant.memory.t_history,1) end;

    # save the SIF and reflectance if there is sunlight
    if _df_dir + _df_dif >= 10
        dfr.BLUE      = MODIS_BLUE(config, spac);
        dfr.EVI       = MODIS_EVI(config, spac);
        dfr.NDVI      = MODIS_NDVI(config, spac);
        dfr.NIR       = MODIS_NIR(config, spac);
        dfr.NIRvI     = MODIS_NIRv(config, spac);
        dfr.NIRvR     = MODIS_NIRvR(config, spac);
        dfr.PAR       = PAR(config, spac);
        dfr.PPAR      = PPAR(spac);
        dfr.RED       = MODIS_RED(config, spac);
        dfr.SIF683    = TROPOMI_SIF683(config, spac);
        dfr.SIF740    = TROPOMI_SIF740(config, spac);
        dfr.SIF757    = OCO2_SIF759(config, spac);
        dfr.SIF771    = OCO2_SIF770(config, spac);
        dfr.ΣSIF      = ΣSIF(spac);
        dfr.ΣSIF_CHL  = ΣSIF_CHL(config, spac);
        dfr.ΣSIF_LEAF = ΣSIF_LEAF(config, spac);
        dfr.ΦD, dfr.ΦF, dfr.ΦN, dfr.ΦP = ΦDFNP(spac);
    else
        dfr.BLUE   = NaN;
        dfr.EVI    = NaN;
        dfr.NDVI   = NaN;
        dfr.NIR    = NaN;
        dfr.NIRvI  = NaN;
        dfr.NIRvR  = NaN;
        dfr.RED    = NaN;
        dfr.ΦD     = NaN;
        dfr.ΦF     = NaN;
        dfr.ΦN     = NaN;
        dfr.ΦP     = NaN;
    end;

    # save water contents and temperatures based on t_on and θ_on
    dfr.MOD_SWC_1 = spac.soils[1].state.θ;
    dfr.MOD_SWC_2 = spac.soils[2].state.θ;
    dfr.MOD_SWC_3 = spac.soils[3].state.θ;
    dfr.MOD_SWC_4 = spac.soils[4].state.θ;

    _tleaf = [leaf.energy.s_aux.t for leaf in spac.plant.leaves];
    dfr.MOD_T_L_MAX  = nanmax(_tleaf);
    dfr.MOD_T_L_MEAN = nanmean(_tleaf);
    dfr.MOD_T_L_MIN  = nanmin(_tleaf);
    dfr.MOD_T_S_1    = spac.soils[1].s_aux.t;
    dfr.MOD_T_S_2    = spac.soils[2].s_aux.t;
    dfr.MOD_T_S_3    = spac.soils[3].s_aux.t;
    dfr.MOD_T_S_4    = spac.soils[4].s_aux.t;

    # save the total flux into the DataFrame
    dfr.BETA  = BETA(spac);
    dfr.F_CO2 = CNPP(spac);
    dfr.F_GPP = GPP(spac);
    dfr.F_H2O = T_VEG(spac);

    return nothing
);

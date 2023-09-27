#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to prescribe parameters from weather drivers
#     2023-Mar-27: prescribe T only if t_on is true, prescribe SWC only is θ_on is true
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: add a controller over rad to make sure it is >= 0
#     2023-Aug-25: make soil and leaf temperature and soil moisture optional
# Bug fixes
#     2023-Aug-26: computed sza in the middle of a time pierod may be > 0 when cumulated radiation is > 0, set it to 88.999
#
#######################################################################################################################################################################################################
"""

    prescribe!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, dfr::DataFrameRow; t_on::Bool = true, θ_on::Bool = true) where {FT}

Prescribe traits and environmental conditions, given
- `config` `SPACConfiguration` type SPAC configuration
- `spac` `MultiLayerSPAC` type SPAC
- `dfr` `DataFrameRow` type weather driver
- `t_on` If true, plant energy budget is on, do not prescribe the temperatures
- `θ_on` If true, soil water budget is on, do not prescribe soil water contents

"""
function prescribe!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, dfr::DataFrameRow; t_on::Bool = true, θ_on::Bool = true) where {FT}
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
    _tleaf = nanmean([_layer.t for _layer in spac.LEAVES]);
    push!(spac.MEMORY.tem, _tleaf);
    if length(spac.MEMORY.tem) > 240 deleteat!(spac.MEMORY.tem,1) end;
    update!(config, spac; t_clm = nanmean(spac.MEMORY.tem));

    # prescribe soil water contents and leaf temperature (for version B1 only)
    if !t_on
        _df_tlf::FT = dfr.T_LEAF;
        _df_ts1::FT = dfr.T_SOIL_1;
        _df_ts2::FT = dfr.T_SOIL_2;
        _df_ts3::FT = dfr.T_SOIL_3;
        _df_ts4::FT = dfr.T_SOIL_4;
        update!(config, spac; t_leaf = max(_df_tar, _df_tlf), t_soils = (_df_ts1, _df_ts2, _df_ts3, _df_ts4));
    end;
    if !θ_on
        _df_sw1::FT = dfr.SWC_1;
        _df_sw2::FT = dfr.SWC_2;
        _df_sw3::FT = dfr.SWC_3;
        _df_sw4::FT = dfr.SWC_4;
        update!(config, spac; swcs = (_df_sw1, _df_sw2, _df_sw3, _df_sw4));
    end;

    # prescribe the precipitation related parameters
    spac.METEO.rain = _df_pcp * ρ_H₂O(FT) / M_H₂O(FT) / 3600;
    spac.METEO.t_precip = _df_tar;

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    _trigger_lai::Bool = !isnan(_df_lai) && (_df_lai != spac.MEMORY.lai);
    _trigger_vcm::Bool = !isnan(_df_vcm) && (_df_vcm != spac.MEMORY.vcm);
    _trigger_chl::Bool = !isnan(_df_chl) && (_df_chl != spac.MEMORY.chl);
    if _trigger_lai
        update!(config, spac; lai = _df_lai, vcmax_expo = 0.3);
        spac.MEMORY.lai = _df_lai;
    end;

    if _trigger_vcm
        update!(config, spac; vcmax = _df_vcm, vcmax_expo = 0.3);
        spac.MEMORY.vcm = _df_vcm;
    end;

    if _trigger_chl
        update!(config, spac; cab = _df_chl, car = _df_chl / 7);
        spac.MEMORY.chl = _df_chl;
    end;

    # update clumping index
    update!(config, spac; ci = _df_cli);

    # update environmental conditions
    for _alayer in spac.AIR
        _alayer.P_AIR = _df_atm;
        update!(_alayer; f_CO₂ = _df_co2, t = _df_tar, vpd = _df_vpd, wind = _df_wnd);
    end;

    # update downward shortwave and longwave radiation
    _in_dir = view(config.SPECTRA.SOLAR_RAD,:,1)' * config.SPECTRA.ΔΛ / 1000;
    _in_dif = view(config.SPECTRA.SOLAR_RAD,:,2)' * config.SPECTRA.ΔΛ / 1000;
    spac.METEO.rad_sw.e_direct  .= view(config.SPECTRA.SOLAR_RAD,:,1) .* max(0,_df_dir) ./ _in_dir;
    spac.METEO.rad_sw.e_diffuse .= view(config.SPECTRA.SOLAR_RAD,:,2) .* max(0,_df_dif) ./ _in_dif;
    spac.METEO.rad_lw = _df_lwr;

    # update solar zenith angle based on the time
    _sza = solar_zenith_angle(spac.LATITUDE, FT(_df_doy));
    spac.ANGLES.sza = (_df_dir + _df_dif > 10) ? min(_sza, 88.999) : _sza;

    return nothing
end


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

    simulation!(wd_tag::String,
                gmdict::Dict{String,Any};
                appending::Bool = false,
                displaying::Bool = false,
                initialial_state::Union{Nothing,Bool} = true,
                saving::Union{Nothing,String} = nothing,
                selection = :)
    simulation!(config::SPACConfiguration{FT},
                spac::MultiLayerSPAC{FT},
                wdf::DataFrame;
                initialial_state::Union{Nothing,Bool} = true,
                saving::Union{Nothing,String} = nothing,
                selection = :) where {FT}

Run simulation on site level, given
- `wd_tag` Weather drive tag such as `wd1`
- `gmdict` GriddingMachine dict for site information
- `appending` If true, append new variables to weather driver when querying the file (set it to true when encountering any errors)
- `displaying` If true, displaying information regarding the steps
- `initialial_state` Initial state of spac: if is a bool, load the first data from the weather driver
- `saving` If is not nothing, save the simulations as a Netcdf file in the working directory; if is nothing, return the simulated result dataframe
- `selection` Run selection of data, default is : (namely 1:end)

The second method can be used to run externally prepared config, spac, and weather driver, given
- `config` SPAC configuration
- `spac` SPAC
- `wdf` Weather driver dataframe

"""
function simulation! end

simulation!(wd_tag::String,
            gmdict::Dict{String,Any};
            appending::Bool = false,
            displaying::Bool = false,
            initialial_state::Union{Nothing,Bool} = true,
            saving::Union{Nothing,String} = nothing,
            selection = :) = (
    _config = spac_config(gmdict);
    _spac = spac(gmdict, _config);
    _wdf = weather_driver(wd_tag, gmdict; appending = appending, displaying = displaying);

    simulation!(_config, _spac, _wdf; initialial_state = initialial_state, saving = saving, selection = selection);

    return isnothing(saving) ? _wdf : nothing
);

simulation!(config::SPACConfiguration{FT},
            spac::MultiLayerSPAC{FT},
            wdf::DataFrame;
            initialial_state::Union{Nothing,Bool} = true,
            saving::Union{Nothing,String} = nothing,
            selection = :) where {FT} = (
    (; DEBUG) = config;

    _wdfr = eachrow(wdf);

    # initialize spac based on initialial_state
    if initialial_state isa Bool
        prescribe!(config, spac, _wdfr[1]; t_on = false, θ_on = false);
    end;

    # iterate through the time steps
    if DEBUG
        for _dfr in _wdfr[selection]
            @show _dfr.ind;
            @time simulation!(config, spac, _dfr);
        end;
    else
        @showprogress for _dfr in _wdfr[selection]
            simulation!(config, spac, _dfr);
        end;
    end;

    # save simulation results to hard drive
    if !isnothing(saving)
        save_nc!(saving, wdf[selection,[_n != "ind" for _n in names(wdf)]]);
    end;

    return nothing
);

simulation!(config::SPACConfiguration{FT},
            spac::MultiLayerSPAC{FT},
            dfr::DataFrameRow;
            n_step::Int = 10,
            δt::Number = 3600
) where {FT} = (
    (; DEBUG, ENABLE_ENERGY_BUDGET, ENABLE_SOIL_WATER_BUDGET) = config;

    # read the data out of dataframe row to reduce memory allocation
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;

    # initialize the integrators
    for _slayer in spac.SOIL.LAYERS
        _slayer.∫∂w∂t_out = 0;
    end;
    for _rlayer in [spac.ROOTS; spac.TRUNK; spac.BRANCHES; spac.LEAVES]
        _rlayer.∫∂w∂t_in = 0;
        _rlayer.∫∂w∂t_out = 0;
    end;

    # prescribe parameters
    prescribe!(config, spac, dfr);

    # run the model
    for _ in 1:n_step
        soil_plant_air_continuum!(config, spac, δt / n_step);
    end;

    # test if the integrated water flow is conserved
    #=
    if DEBUG
        _sum_soil_out = sum([_slayer.∫∂w∂t_out for _slayer in spac.SOIL.LAYERS]);
        _sum_root_in = sum([_rlayer.∫∂w∂t_in for _rlayer in spac.ROOTS]) / spac.SOIL.AREA;
        _sum_root_out = sum([_rlayer.∫∂w∂t_out for _rlayer in spac.ROOTS]) / spac.SOIL.AREA;
        _sum_trnk_in = spac.TRUNK.∫∂w∂t_in / spac.SOIL.AREA;
        _sum_trnk_out = spac.TRUNK.∫∂w∂t_out / spac.SOIL.AREA;
        _sum_stem_in = sum([_slayer.∫∂w∂t_in for _slayer in spac.BRANCHES]) / spac.SOIL.AREA;
        _sum_stem_out = sum([_slayer.∫∂w∂t_out for _slayer in spac.BRANCHES]) / spac.SOIL.AREA;
        _sum_leaf_in = [_slayer.∫∂w∂t_in for _slayer in spac.LEAVES]' * spac.CANOPY.δlai;
        _sum_leaf_out = [_clayer.∫∂w∂t_out for _clayer in spac.LEAVES]' * spac.CANOPY.δlai;
        @info "Debugging" _sum_soil_out _sum_root_in _sum_root_out _sum_trnk_in _sum_trnk_out _sum_stem_in _sum_stem_out _sum_leaf_in _sum_leaf_out spac.METEO.rain;
    end;
    =#

    # save the SIF and reflectance if there is sunlight
    if _df_dir + _df_dif >= 10
        dfr.BLUE      = MODIS_BLUE(config, spac);
        dfr.EVI       = MODIS_EVI(config, spac);
        dfr.NDVI      = MODIS_NDVI(config, spac);
        dfr.NIR       = MODIS_NIR(config, spac);
        dfr.NIRvI     = MODIS_NIRv(config, spac);
        dfr.NIRvR     = MODIS_NIRvR(config, spac);
        dfr.PAR       = spac.CANOPY.RADIATION.par_in;
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

        # display the debug information
        if DEBUG
            if any(isnan, (dfr.BLUE, dfr.EVI, dfr.NDVI, dfr.NIR, dfr.NIRvI, dfr.NIRvR, dfr.PAR, dfr.PPAR, dfr.RED, dfr.SIF683, dfr.SIF740, dfr.SIF757, dfr.SIF771, dfr.ΦD, dfr.ΦF, dfr.ΦN, dfr.ΦP))
                @info "Debugging" dfr.BLUE dfr.EVI dfr.NDVI dfr.NIR dfr.NIRvI dfr.NIRvR dfr.PAR dfr.PPAR dfr.RED dfr.SIF683 dfr.SIF740 dfr.SIF757 dfr.SIF771 dfr.ΦD dfr.ΦF dfr.ΦN dfr.ΦP;
                error("NaN detected when computing remote sensing variables");
            end;
        end;
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
    if ENABLE_SOIL_WATER_BUDGET
        dfr.MOD_SWC_1 = spac.SOIL.LAYERS[1].θ;
        dfr.MOD_SWC_2 = spac.SOIL.LAYERS[2].θ;
        dfr.MOD_SWC_3 = spac.SOIL.LAYERS[3].θ;
        dfr.MOD_SWC_4 = spac.SOIL.LAYERS[4].θ;

        if DEBUG
            if any(isnan, (dfr.MOD_SWC_1, dfr.MOD_SWC_2, dfr.MOD_SWC_3, dfr.MOD_SWC_4))
                @info "Debugging" dfr.MOD_SWC_1 dfr.MOD_SWC_2 dfr.MOD_SWC_3 dfr.MOD_SWC_4;
                error("NaN detected when computing soil water contents");
            end;
        end;
    end;
    if ENABLE_ENERGY_BUDGET
        _tleaf = [_leaf.t for _leaf in spac.LEAVES];
        dfr.MOD_T_L_MAX  = nanmax(_tleaf);
        dfr.MOD_T_L_MEAN = nanmean(_tleaf);
        dfr.MOD_T_L_MIN  = nanmin(_tleaf);
        dfr.MOD_T_S_1    = spac.SOIL.LAYERS[1].t;
        dfr.MOD_T_S_2    = spac.SOIL.LAYERS[2].t;
        dfr.MOD_T_S_3    = spac.SOIL.LAYERS[3].t;
        dfr.MOD_T_S_4    = spac.SOIL.LAYERS[4].t;

        if DEBUG
            if any(isnan, (dfr.MOD_T_L_MAX, dfr.MOD_T_L_MEAN, dfr.MOD_T_L_MIN, dfr.MOD_T_S_1, dfr.MOD_T_S_2, dfr.MOD_T_S_3, dfr.MOD_T_S_4))
                @info "Debugging" dfr.MOD_T_L_MAX dfr.MOD_T_L_MEAN dfr.MOD_T_L_MIN dfr.MOD_T_S_1 dfr.MOD_T_S_2 dfr.MOD_T_S_3 dfr.MOD_T_S_4;
                error("NaN detected when computing soil and leaf temperatures");
            end;
        end;
    end;

    # save the total flux into the DataFrame
    dfr.BETA  = BETA(spac);
    dfr.F_CO2 = CNPP(spac);
    dfr.F_GPP = GPP(spac);
    dfr.F_H2O = T_VEG(spac);

    if DEBUG
        if any(isnan, (dfr.F_CO2, dfr.F_GPP, dfr.F_H2O))
            @info "Debugging" dfr.BETA dfr.F_CO2 dfr.F_GPP dfr.F_H2O;
            error("NaN detected when computing fluxes");
        end;
    end;

    return nothing
);

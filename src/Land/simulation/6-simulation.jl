"""

    simulation!(settings::Union{Dict,OrderedDict}, lat::Number, lon::Number, year::Int; saving::Union{Nothing,String} = nothing)
    simulation!(settings::Union{Dict,OrderedDict}, gmd::Dict{String,Any}; saving::Union{Nothing,String} = nothing)

Run simulation on site level, given
- `settings` Dictionary of settings
- `lat` Latitude of the site
- `lon` Longitude of the site
- `year` Year of the simulation
- `saving` If is not nothing, save the simulations as a Netcdf file in the working directory; if is nothing, return the simulated result dataframe
- `gmd` GriddingMachine dict for site information

"""
function simulation! end;

simulation!(settings::Union{Dict,OrderedDict}, lat::Number, lon::Number, year::Int; saving::Union{Nothing,String} = nothing) =
    simulation!(settings, grid_dict(LandDatasetLabels(settings["GM_VERSION"], year), lat, lon); saving = saving);

simulation!(settings::Union{Dict,OrderedDict}, gmd::Dict{String,Any}; saving::Union{Nothing,String} = nothing) = (
    wd = grid_weather(WeatherDriverLabels(settings["WD_VERSION"], gmd["YEAR"]), gmd["LATITUDE"], gmd["LONGITUDE"]);
    sd = parameters_to_save(settings["VARIABLES_TO_SAVE"]);
    config = site_config(settings);
    spac = site_spac(config, gmd);
    driver = site_driver_tuple(gmd, wd);
    results = site_result_tuple(spac, wd, sd);

    return simulation!(config, spac, driver, results; saving = saving, saving_setting = sd, selection = settings["SIMULATION_PERIOD"], δt = settings["TIME_STEP"]);
);

simulation!(config::SPACConfig{FT},
            spac::BulkSPAC{FT},
            driver::NamedTuple,
            results::NamedTuple;
            saving::Union{Nothing,String} = nothing,
            saving_setting::Vector{ParameterFunctionMapper} = parameters_to_save(),
            selection = :,
            δt::Number = 3600) where {FT} = (
    (; MESSAGE_LEVEL) = config.CONFIG_INFO;

    # initialize spac based on initialize_state for the first time step
    prescribe!(config, spac, driver, 1; initialize_state = true);

    # iterate through the time steps
    if MESSAGE_LEVEL == 0
        for idx in eachindex(driver.FDOY)[selection]
            simulate_step!(config, spac, driver, results, idx, saving_setting, δt);
        end;
    elseif MESSAGE_LEVEL == 1
        @showprogress for idx in eachindex(driver.FDOY)[selection]
            simulate_step!(config, spac, driver, results, idx, saving_setting, δt);
        end;
    elseif MESSAGE_LEVEL == 2
        for idx in eachindex(driver.FDOY)[selection]
            print("\rRunning simulation for $(lpad(idx,4," ")) out of $(lpad(length(driver.FDOY),4," "))...");
            simulate_step!(config, spac, driver, results, idx, saving_setting, δt);
        end;
    else
        error("MESSAGE_LEVEL should be 0, 1, or 2");
    end;

    # save simulation results to hard drive
    df = DataFrame(results);

    return isnothing(saving) ? df[selection, names(df)] : save_nc!(saving, df[selection, names(df)])
);


# function to simulate a single time step
function simulate_step!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, driver::NamedTuple, results::NamedTuple, ind::Int, saving_setting::Vector{ParameterFunctionMapper}, δt::Number) where {FT}
    # prescribe parameters
    prescribe!(config, spac, driver, ind);

    # run the model
    soil_plant_air_continuum!(config, spac, δt);
    push_t_history!(config, spac);

    # save the results
    for lpm in saving_setting
        if lpm.to_save
            results[Symbol(lpm.name)][ind] = lpm.func(config, spac, lpm.params...);
        end;
    end;

    return nothing
end;

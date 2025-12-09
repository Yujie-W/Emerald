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
#     2024-Mar-07: add fields for saved parameters and simulations in the dataframe (as grid_weather_driver did not do it)
#     2024-Aug-05: use saving_dict to determine which variables to save
#     2024-Aug-05: save plant hydraulics health status
#     2024-Aug-05: add method to use externally prepared config, spac, and weather driver (will process the dataframe to NamedTuple)
#     2024-Aug-05: add option to save soil water potential
#     2025-Sep-04: remove the option to prescribe leaf temperature (need to spin up in the future)
#
#######################################################################################################################################################################################################
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

    return simulation!(config, spac, driver, results; saving = saving, saving_dict = sd, selection = settings["SIMULATION_PERIOD"], δt = settings["TIME_STEP"]);
);

simulation!(config::SPACConfig{FT},
            spac::BulkSPAC{FT},
            driver::NamedTuple,
            results::NamedTuple;
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Bool} = parameters_to_save(),
            selection = :,
            δt::Number = 3600) where {FT} = (
    (; MESSAGE_LEVEL) = config.CONFIG_INFO;

    # initialize spac based on initialize_state for the first time step
    prescribe!(config, spac, driver, 1; initialize_state = true);

    # iterate through the time steps
    if MESSAGE_LEVEL == 0
        for idx in eachindex(driver.FDOY)[selection]
            simulation!(config, spac, driver, results, idx; saving_dict = saving_dict, δt = δt);
        end;
    elseif MESSAGE_LEVEL == 1
        @showprogress for idx in eachindex(driver.FDOY)[selection]
            simulation!(config, spac, driver, results, idx; saving_dict = saving_dict, δt = δt);
        end;
    elseif MESSAGE_LEVEL == 2
        for idx in eachindex(driver.FDOY)[selection]
            print("\rRunning simulation for $(lpad(idx,4," ")) out of $(lpad(length(driver.FDOY),4," "))...");
            simulation!(config, spac, driver, results, idx; saving_dict = saving_dict, δt = δt);
        end;
    else
        error("MESSAGE_LEVEL should be 0, 1, or 2");
    end;

    # save simulation results to hard drive
    df = DataFrame(results);

    return isnothing(saving) ? df[selection, names(df)] : save_nc!(saving, df[selection, names(df)])
);

simulation!(config::SPACConfig{FT},
            spac::BulkSPAC{FT},
            driver::NamedTuple,
            results::NamedTuple,
            ind::Int;
            saving_dict::Dict{String,Bool} = parameters_to_save(),
            δt::Number = 3600) where {FT} = (
    # prescribe parameters
    prescribe!(config, spac, driver, ind);

    # run the model
    soil_plant_air_continuum!(config, spac, δt);
    push_t_history!(config, spac);

    # save the results
    save_fields!(config, spac, results, ind; saving_dict = saving_dict);

    return nothing
);

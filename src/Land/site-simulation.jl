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

    simulation!(wd_tag::String,
                gm_dict::Dict{String,Any};
                appending::Bool = false,
                saving::Union{Nothing,String} = nothing,
                saving_dict::Dict{String,Bool} = SAVING_DICT,
                selection = :)
    simulation!(config::SPACConfig{FT},
                spac::BulkSPAC{FT},
                df::DataFrame;
                saving::Union{Nothing,String} = nothing,
                saving_dict::Dict{String,Bool} = SAVING_DICT,
                selection = :) where {FT}

Run simulation on site level, given
- `wd_tag` Weather drive tag such as `wd1`
- `gm_dict` GriddingMachine dict for site information
- `appending` If true, append new variables to weather driver when querying the file (set it to true when encountering any errors)
- `saving` If is not nothing, save the simulations as a Netcdf file in the working directory; if is nothing, return the simulated result dataframe
- `selection` Run selection of data, default is : (namely 1:end;)

The second method can be used to run externally prepared config, spac, and weather driver, given
- `config` SPAC configuration
- `spac` SPAC
- `df` Weather driver dataframe

"""
function simulation! end;

simulation!(gm_tag::String,
            wd_tag::String,
            lat::Number,
            lon::Number,
            year::Int;
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Bool} = DEFAULT_SAVING_DICT,
            selection = :) = simulation!(grid_dict(LandDatasetLabels(gm_tag, year), lat, lon), wd_tag, lat, lon, year; saving = saving, saving_dict = saving_dict, selection = selection);

simulation!(gmd::Dict{String,Any},
            wd_tag::String,
            lat::Number,
            lon::Number,
            year::Int;
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Bool} = DEFAULT_SAVING_DICT,
            selection = :) = (
    wd = grid_weather(WeatherDriverLabels(wd_tag, year), lat, lon);
    config = site_config(gmd);
    spac = site_spac(config, gmd);
    driver = site_driver_tuple(gmd, wd);
    results = site_result_tuple(spac, wd, parameters_to_save());

    return simulation!(config, spac, driver, results; saving = saving, saving_dict = saving_dict, selection = selection);
);

simulation!(config::SPACConfig{FT},
            spac::BulkSPAC{FT},
            driver::NamedTuple,
            results::NamedTuple;
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Bool} = DEFAULT_SAVING_DICT,
            selection = :) where {FT} = (
    (; MESSAGE_LEVEL) = config.CONFIG_INFO;

    # initialize spac based on initialize_state for the first time step
    prescribe!(config, spac, driver, 1; initialize_state = true);

    # iterate through the time steps
    if MESSAGE_LEVEL == 0
        for idx in eachindex(driver.FDOY)[selection]
            simulation!(config, spac, driver, results, idx; saving_dict = saving_dict);
        end;
    elseif MESSAGE_LEVEL == 1
        @showprogress for idx in eachindex(driver.FDOY)[selection]
            simulation!(config, spac, driver, results, idx; saving_dict = saving_dict);
        end;
    elseif MESSAGE_LEVEL == 2
        for idx in eachindex(driver.FDOY)[selection]
            print("Running simulation for $idx out of $(length(driver.FDOY))...");
            simulation!(config, spac, driver, results, idx; saving_dict = saving_dict);
        end;
    else
        error("MESSAGE_LEVEL should be 0, 1, or 2");
    end;

    # save simulation results to hard drive
    if !isnothing(saving)
        df = DataFrame(results);
        save_nc!(saving, df[selection, names(df)]);
    end;

    return nothing
);

simulation!(config::SPACConfig{FT},
            spac::BulkSPAC{FT},
            driver::NamedTuple,
            results::NamedTuple,
            ind::Int;
            saving_dict::Dict{String,Bool} = DEFAULT_SAVING_DICT,
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

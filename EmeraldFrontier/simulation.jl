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
#
#######################################################################################################################################################################################################
"""

    simulation!(wd_tag::String,
                gm_dict::Dict{String,Any};
                appending::Bool = false,
                initialize_state::Union{Nothing,Bool} = true,
                saving::Union{Nothing,String} = nothing,
                saving_dict::Dict{String,Any} = SAVING_DICT,
                selection = :)
    simulation!(config::SPACConfiguration{FT},
                spac::BulkSPAC{FT},
                df::DataFrame;
                initialize_state::Union{Nothing,Bool} = true,
                saving::Union{Nothing,String} = nothing,
                saving_dict::Dict{String,Any} = SAVING_DICT,
                selection = :) where {FT}

Run simulation on site level, given
- `wd_tag` Weather drive tag such as `wd1`
- `gm_dict` GriddingMachine dict for site information
- `appending` If true, append new variables to weather driver when querying the file (set it to true when encountering any errors)
- `initialize_state` Initial state of spac: if is a bool, load the first data from the weather driver
- `saving` If is not nothing, save the simulations as a Netcdf file in the working directory; if is nothing, return the simulated result dataframe
- `selection` Run selection of data, default is : (namely 1:end;)

The second method can be used to run externally prepared config, spac, and weather driver, given
- `config` SPAC configuration
- `spac` SPAC
- `df` Weather driver dataframe

"""
function simulation! end;

simulation!(wd_tag::String,
            gm_dict::Dict{String,Any};
            appending::Bool = false,
            initialize_state::Union{Nothing,Bool} = true,
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Any} = SAVING_DICT,
            selection = :) = (
    config = spac_config(gm_dict);
    spac = grid_spac(config, gm_dict);
    df = grid_weather_driver(wd_tag, gm_dict; appending = appending);

    return simulation!(config, spac, df; initialize_state = initialize_state, saving = saving, saving_dict = saving_dict, selection = selection);
);

simulation!(config::SPACConfiguration{FT},
            spac::BulkSPAC{FT},
            df::DataFrame;
            initialize_state::Union{Nothing,Bool} = true,
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Any} = SAVING_DICT,
            selection = :) where {FT} = (
    # convert the DataFrame to NamedTuple with new fields
    wdf = prepare_wdf(spac, df; saving_dict = saving_dict);

    simulation!(config, spac, wdf; initialize_state = initialize_state, saving = saving, saving_dict = saving_dict, selection = selection);

    return isnothing(saving) ? DataFrame(wdf) : nothing
);

simulation!(config::SPACConfiguration{FT},
            spac::BulkSPAC{FT},
            wdf::NamedTuple;
            initialize_state::Union{Nothing,Bool} = true,
            saving::Union{Nothing,String} = nothing,
            saving_dict::Dict{String,Any} = SAVING_DICT,
            selection = :) where {FT} = (
    (; MESSAGE_LEVEL) = config;

    # initialize spac based on initialize_state
    prescribe!(config, spac, wdf, 1; initialize_state = initialize_state);

    # iterate through the time steps
    if MESSAGE_LEVEL == 0
        for idx in eachindex(wdf.FDOY)[selection]
            simulation!(config, spac, wdf, idx; saving_dict = saving_dict);
        end;
    elseif MESSAGE_LEVEL == 1
        @showprogress for idx in eachindex(wdf.FDOY)[selection]
            simulation!(config, spac, wdf, idx; saving_dict = saving_dict);
        end;
    elseif MESSAGE_LEVEL == 2
        for idx in eachindex(wdf.FDOY)[selection]
            @show wdf.ind[idx];
            simulation!(config, spac, wdf, idx; saving_dict = saving_dict);
        end;
    else
        error("MESSAGE_LEVEL should be 0, 1, or 2");
    end;

    # save simulation results to hard drive
    if !isnothing(saving)
        df = DataFrame(wdf);
        save_nc!(saving, df[selection, [n != "ind" for n in names(df)]]);
    end;

    return nothing
);

simulation!(config::SPACConfiguration{FT},
            spac::BulkSPAC{FT},
            wdf::NamedTuple,
            ind::Int;
            saving_dict::Dict{String,Any} = SAVING_DICT,
            δt::Number = 3600) where {FT} = (
    # prescribe parameters
    prescribe!(config, spac, wdf, ind);

    # run the model
    soil_plant_air_continuum!(config, spac, δt);
    push_t_history!(config, spac);

    # save the results
    save_fields!(config, spac, wdf, ind; saving_dict = saving_dict);

    return nothing
);

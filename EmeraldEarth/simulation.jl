#######################################################################################################################################################################################################
#
# Changes to this file
# General
#     2024-Feb-28: redesign the global simulations
#     2024-Feb-28: add option single_thread for debugging purpose
#
#######################################################################################################################################################################################################
"""

    global_simulations!(
                gm_mat::Matrix{Union{Nothing,Dict{String,Any}}},
                wd_mat::Matrix{Dict{String,Any}},
                st_mat::Matrix{Union{Nothing,BulkSPACStates{FT}}};
                single_thread::Bool = false) where {FT}

Run global simulations per time step, given
- `gm_mat` GriddingMachine inputs matrix
- `wd_mat` Weather drivers matrix
- `st_mat` Initial states matrix
- `single_thread` Run the simulations in single thread mode for debugging purpose

"""
function global_simulations!(
            gm_mat::Matrix{Union{Nothing,Dict{String,Any}}},
            wd_mat::Matrix{Dict{String,Any}},
            st_mat::Matrix{Union{Nothing,BulkSPACStates{FT}}};
            single_thread::Bool = false) where {FT}
    if single_thread
        for i in axes(gm_mat, 1), j in axes(gm_mat, 2)
            @tinfo "Running simulation for grid $(i), $(j)";
            gm_dict = gm_mat[i,j];
            wd_dict = wd_mat[i,j];
            state = st_mat[i,j];
            new_state = grid_simulation!(gm_dict, wd_dict, state);
        end;

        return nothing
    end;

    # run the simulations in multiple threads
    states = @showprogress pmap(grid_simulation!, gm_mat, wd_mat, st_mat);

    return states
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#     2023-Mar-13: add state mat as an input
#     2023-Apr-13: add CACHE_CONFIG in function
#     2023-Jun-15: add some code for debugging (not used for actual simulations)
#     2023-Jun-15: add option to display process information
#     2024-Feb-29: if wd_dict contains NaN, return nothing (hereafter)
#
#######################################################################################################################################################################################################
"""

    grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::BulkSPACStates)
    grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::Nothing)
    grid_simulation!(gm_dict::Nothing, wd_dict::Dict{String,Any}, state::Nothing)

Run simulations on SPAC, given
- `gm_dict` GriddingMachine inputs
- `wd_dict` Weather drivers
- `state` Initial states

"""
function grid_simulation! end;

grid_simulation!(gm_dict::Nothing, wd_dict::Dict{String,Any}, state::Nothing) = nothing;

grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::Nothing) = nothing;

grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::BulkSPACStates) = (
    # if wd_dict contains NaN, return nothing
    if dict_contains_nan(wd_dict)
        #@warn "Weather driver contains NaN, skipping the simulation..." gm_dict["LATITUDE"] gm_dict["LONGITUDE"];
        return nothing
    end;

    # continue only if wd_dict does not have NaN
    spac = grid_spac(CACHE_CONFIG, gm_dict);
    prescribe_gm_wd_data!(CACHE_CONFIG, spac, gm_dict, wd_dict);
    initialize_spac!(CACHE_CONFIG, spac, state);
    soil_plant_air_continuum!(CACHE_CONFIG, spac, 3600);

    # update the temperature history
    mean_tleaf = nanmean([l.energy.s_aux.t for l in spac.plant.leaves]);
    push!(spac.plant.memory.t_history, mean_tleaf);
    if length(spac.plant.memory.t_history) > 240 deleteat!(spac.plant.memory.t_history,1) end;

    return BulkSPACStates(spac)
);

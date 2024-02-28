function global_simulation! end;

global_simulation!(mat_gm::Matrix{Union{Nothing,Dict{String,Any}}}, mat_wd::Matrix{Union{Nothing,Dict{String,Any}}}, mat_st::Matrix{Union{Nothing,BulkSPACStates{FT}}} where FT) = (
    states = @showprogress pmap(grid_simulation!, mat_gm, mat_wd, mat_st);

    return states
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#     2023-Mar-13: add state mat as an input
#     2023-Apr-13: add CACHE_CONFIG in function
#     2023-Jun-15: add some code for debugging (not used for actual simulations)
#     2023-Jun-15: add option to display process information
#
#######################################################################################################################################################################################################
"""

    grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::BulkSPACStates)
    grid_simulation!(gm_dict::Nothing, wd_dict::Dict{String,Any}, state::Nothing)

Run simulations on SPAC, given
- `gm_dict` GriddingMachine inputs
- `wd_dict` Weather drivers
- `state` Initial states

"""
function grid_simulation! end;

grid_simulation!(gm_dict::Nothing, wd_dict::Dict{String,Any}, state::Nothing) = nothing;

grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::BulkSPACStates) = (
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

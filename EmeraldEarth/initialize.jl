#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-23: add function to initialize the state based on the gridding machine, weather driver, and initial states matrices
#
#######################################################################################################################################################################################################
"""

    initial_states(mat_gm::Matrix{Union{Nothing,Dict{String,Any}}}, mat_wd::Matrix{Dict{String,Any}}, mat_ss::Matrix{Dict{String,Any}})

Initialize the state using multiple threads, given
- `mat_gm` Matrix of GriddingMachine inputs
- `mat_wd` Matrix of weather drivers
- `mat_ss` Matrix of initial states

"""
function initial_states(mat_gm::Matrix{Union{Nothing,Dict{String,Any}}}, mat_wd::Matrix{Dict{String,Any}}, mat_ss::Matrix{Dict{String,Any}})
    states = @showprogress pmap(initial_state, mat_gm, mat_wd, mat_ss);

    return states
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-23: add function to initialize the state based on the gridding machine, weather driver, and initial state
#     2024-Feb-23: set up SAI as well
#     2024-Feb-28: move some prescribe gm, wd, and ss functions to EmeraldData
#
#######################################################################################################################################################################################################
"""

    initial_state(gm_dict::Nothing, wd_dict::Dict{String,Any}, ss_dict::Dict{String,Any})
    initial_state(gm_dict::Union{Nothing,Dict{String,Any}}, wd_dict::Dict{String,Any}, ss_dict::Dict{String,Any})

Initialize the state, given
- `gm_dict` GriddingMachine inputs
- `wd_dict` Weather drivers
- `ss_dict` Initial states

"""
function initial_state end;

initial_state(gm_dict::Nothing, wd_dict::Dict{String,Any}, ss_dict::Dict{String,Any}) = nothing;

initial_state(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, ss_dict::Dict{String,Any}) = (
    FT = gm_dict["FT"];
    spac = grid_spac(CACHE_CONFIG, gm_dict);
    prescribe_gm_wd_data!(CACHE_CONFIG, spac, gm_dict, wd_dict, ss_dict);

    # initialize the spac with non-saturated soil
    initialize_spac!(CACHE_CONFIG, spac);
    soil_plant_air_continuum!(CACHE_CONFIG, spac, FT(0));

    return BulkSPACStates(spac)
);

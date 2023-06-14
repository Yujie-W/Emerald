#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#     2023-Mar-13: add state mat as an input
#     2023-Apr-13: add CACHE_CONFIG in function
#
#######################################################################################################################################################################################################
"""

    simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}}, wd_mat::Matrix{Union{Nothing,Dict{String,Any}}}, state_mat::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}) where {FT<:AbstractFloat}

Run simulations on SPAC, given
- `gm_mat` Matrix of GriddingMachine inputs
- `wd_mat` Matrix of weather drivers
- `state_mat` Matrix of state variable struct

"""
function simulation! end

simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}}, wd_mat::Matrix{Union{Nothing,Dict{String,Any}}}, state_mat::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}) where {FT<:AbstractFloat} = (
    @tinfo "Running the global simulations in multiple threads...";
    _states = @showprogress pmap(simulation!, gm_mat, wd_mat, state_mat);

    return _states
);

simulation!(gm_params::Nothing, wd_params::Nothing, state::Nothing) = nothing;

simulation!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MultiLayerSPACState{FT}}) where {FT<:AbstractFloat} = (
    synchronize_cache!(gm_params, wd_params, state);
    for _ in 1:10
        soil_plant_air_continuum!(CACHE_SPAC, CACHE_CONFIG, 360; p_on = false, t_on = false, θ_on = false);
    end;
    spac_state!(CACHE_SPAC, CACHE_STATE);

    # if wd_params["RAD_DIR"] > 500
    #     @info "Debugging" CACHE_SPAC.LATITUDE CACHE_SPAC.LONGITUDE GPP(CACHE_SPAC) PPAR(CACHE_SPAC);
    # end;

    return CACHE_STATE
);

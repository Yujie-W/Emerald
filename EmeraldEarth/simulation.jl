#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#     2023-Mar-13: add state mat as an input
#     2023-Apr-13: add CACHE_CONFIG in function
#     2023-Jun-15: add some code for debugging (not used for actual simulations)
#     2023-Jun-15: add option to display process information
# To do
#     TODO: use integrated value for GPP, ET, etc
#
#######################################################################################################################################################################################################
"""

    simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}},
                wd_mat::Matrix{Union{Nothing,Dict{String,Any}}},
                state_mat::Matrix{Union{Nothing}};
                displaying::Bool = true
    )

Run simulations on SPAC, given
- `gm_mat` Matrix of GriddingMachine inputs
- `wd_mat` Matrix of weather drivers
- `state_mat` Matrix of state variable struct
- `displaying` Whether to display information regarding process

"""
function simulation! end;

simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}},
            wd_mat::Matrix{Union{Nothing,Dict{String,Any}}},
            state_mat::Matrix{Union{Nothing}};
            displaying::Bool = true
) = (
#    @tinfo "Debugging the code using one core...";
#    @showprogress for i in eachindex(gm_mat)
#        simulation!(gm_mat[i], wd_mat[i], state_mat[i]);
#        if !(0 <= CACHE_STATE.beta <= 1) ||
#           !(CACHE_STATE.csif >= 0) ||
#           !(CACHE_STATE.etr >= 0) ||
#           !(CACHE_STATE.gpp >= 0) ||
#           !(CACHE_STATE.oco_sif₇₅₉ >= 0) ||
#           !(CACHE_STATE.oco_sif₇₇₀ >= 0) ||
#           !(CACHE_STATE.par >= 0) ||
#           !(CACHE_STATE.ppar >= 0) ||
#           !(CACHE_STATE.tropomi_sif₆₈₃ >= 0) ||
#           !(CACHE_STATE.tropomi_sif₇₄₀ >= 0) ||
#           isnan(CACHE_STATE.transpiration)
#            @warn "An error occurred 1!"; break;
#        end;
#        if CACHE_SPAC.canopy.sun_geometry.state.sza < 89
#            if isnan(CACHE_STATE.modis_evi) ||
#               isnan(CACHE_STATE.modis_ndvi) ||
#               isnan(CACHE_STATE.modis_nirv)
#                @warn "An error occurred 2!"; break;
#            end;
#        else
#            if !isnan(CACHE_STATE.modis_evi) ||
#               !isnan(CACHE_STATE.modis_ndvi) ||
#               !isnan(CACHE_STATE.modis_nirv)
#                @warn "An error occurred 3!"; break;
#            end;
#        end;
#    end;
#
#    return state_mat

    if displaying
        @tinfo "Running the global simulations in multiple threads...";
        _states = @showprogress pmap(simulation!, gm_mat, wd_mat, state_mat);
    else
        _states = pmap(simulation!, gm_mat, wd_mat, state_mat);
    end;

    return _states
);

simulation!(gm_params::Nothing, wd_params::Nothing, state::Nothing) = nothing;

simulation!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing}) = (
    synchronize_cache!(gm_params, wd_params, state);
    soil_plant_air_continuum!(CACHE_CONFIG, CACHE_SPAC, 3600);
    #spac_state!(CACHE_CONFIG, CACHE_SPAC, CACHE_STATE);

    return CACHE_STATE
);

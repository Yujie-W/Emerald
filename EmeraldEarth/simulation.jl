#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: add step to synchronize state variables into CACHE_SPAC
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-15: make sure prescribed soil parameters are not NaN and rad is >= 0
#
#######################################################################################################################################################################################################
"""

    prescribe_traits_environment!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any})

Prescribe traits and environmental conditions to SPAC, given
- `gm_dict` Dict for GriddingMachine parameters
- `wd_dict` Dict for weather drivers

"""
function prescribe_traits_environment!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any})
    FT = gm_dict["FT"];

    # update environmental conditions
    for air in CACHE_SPAC.airs
        air.state.p_air = wd_dict["P_ATM"];
        prescribe_air!(air; f_CO₂ = gm_dict["CO2"], t = wd_dict["T_AIR"], vpd = wd_dict["VPD"], wind = wd_dict["WIND"]);
    end;

    # update shortwave and longwave radiation
    ref_dir = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1)'  * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    ref_dif = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2)' * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    CACHE_SPAC.meteo.rad_sw.e_dir .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1) .* max(0,wd_dict["RAD_DIR"]) ./ ref_dir;
    CACHE_SPAC.meteo.rad_sw.e_dif .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2) .* max(0,wd_dict["RAD_DIF"]) ./ ref_dif;
    CACHE_SPAC.meteo.rad_lw = wd_dict["RAD_LW"];
    sza = solar_zenith_angle(CACHE_SPAC.info.lat, FT(wd_dict["FDOY"]));
    CACHE_SPAC.canopy.sun_geometry.state.sza = (wd_dict["RAD_DIR"] + wd_dict["RAD_DIF"] > 10) ? min(sza, 88.999) : sza;

    # update t_clm to make Vcmax25 and Jmax25 TD temperature dependent
    prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; t_clm = nanmean(CACHE_SPAC.plant.memory.t_history));

    # synchronize LAI, CHL, and CI
    # TODO: use memory to double check if these need to be updated
    iday = Int(floor(wd_dict["INDEX"] / 24)) + 1;
    chl = query_griddingmachine_data(gm_dict["CHLOROPHYLL"], gm_dict["YEAR"], iday);
    ci = query_griddingmachine_data(gm_dict["CLUMPING"], gm_dict["YEAR"], iday);
    lai = query_griddingmachine_data(gm_dict["LAI"], gm_dict["YEAR"], iday);
    vcm = query_griddingmachine_data(gm_dict["VCMAX25"], gm_dict["YEAR"], iday);
    prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; ci = ci, lai = lai, vcmax = vcm, vcmax_expo = 0.3);
    prescribe_traits!(CACHE_CONFIG, CACHE_SPAC; cab = chl, car = chl / 7);

    # initialize the spac with non-saturated soil
    # TODO: do something here
    # initialize_spac!(CACHE_CONFIG, CACHE_SPAC);

    return nothing
end;



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
function grid_simulation! end;

grid_simulation!(gm_dict::Nothing, wd_dict::Dict{String,Any}, state::Nothing) = nothing;

grid_simulation!(gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, state::BulkSPACStates) = (
    initialize_spac!(CACHE_CONFIG, CACHE_SPAC, state);
    prescribe_traits_environment!(gm_dict, wd_dict);
    soil_plant_air_continuum!(CACHE_CONFIG, CACHE_SPAC, 3600);
    mean_tleaf = nanmean([l.energy.s_aux.t for l in CACHE_SPAC.plant.leaves]);
    push!(CACHE_SPAC.plant.memory.t_history, mean_tleaf);
    if length(CACHE_SPAC.plant.memory.t_history) > 240 deleteat!(CACHE_SPAC.plant.memory.t_history,1) end;

    sync_state!(CACHE_SPAC, state);

    return state
);

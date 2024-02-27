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
    #
    # TODO: add support to C4 photosynthesis
    #
    FT = gm_dict["FT"];
    zc = max(FT(0.1), gm_dict["CANOPY_HEIGHT"]);
    spac = BulkSPAC(
                CACHE_CONFIG;
                air_bounds = collect(0:21) * zc / 20,
                elevation = gm_dict["ELEVATION"],
                latitude = gm_dict["LATITUDE"],
                longitude = gm_dict["LONGITUDE"],
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                plant_zs = [-2, zc/2, zc]);
    spac.soil_bulk.trait.color = gm_dict["SOIL_COLOR"];
    for i in eachindex(spac.plant.leaves)
        spac.plant.leaves[i].bio.trait.lma = gm_dict["LMA"];
        spac.plant.leaves[i].flux.trait.stomatal_model = deepcopy(CACHE_SPAC.plant.leaves[i].flux.trait.stomatal_model);
        spac.plant.leaves[i].flux.trait.stomatal_model.G1 = gm_dict["G1_MEDLYN_C3"];
    end;

    # set up SAI
    prescribe_traits!(CACHE_CONFIG, spac; sai = gm_dict["SAI"]);

    # update soil type information per layer
    for i in eachindex(spac.soils)
        # TODO: add a line to parameterize K_MAX
        # TODO: fix these later with better data source
        if !isnan(gm_dict["SOIL_α"][i]) && !isnan(gm_dict["SOIL_N"][i]) && !isnan(gm_dict["SOIL_ΘR"][i]) && !isnan(gm_dict["SOIL_ΘS"][i])
            spac.soils[i].trait.vc.α = gm_dict["SOIL_α"][i];
            spac.soils[i].trait.vc.N = gm_dict["SOIL_N"][i];
            spac.soils[i].trait.vc.M = 1 - 1 / spac.soils[i].trait.vc.N;
            spac.soils[i].trait.vc.Θ_RES = gm_dict["SOIL_ΘR"][i];
            spac.soils[i].trait.vc.Θ_SAT = gm_dict["SOIL_ΘS"][i];
        end;
    end;

    # update environmental conditions
    for air in spac.airs
        air.state.p_air = wd_dict["P_ATM"];
        prescribe_air!(air; f_CO₂ = gm_dict["CO2"], t = wd_dict["T_AIR"], vpd = wd_dict["VPD"], wind = wd_dict["WIND"]);
    end;

    # update shortwave and longwave radiation
    ref_dir = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1)'  * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    ref_dif = view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2)' * CACHE_CONFIG.SPECTRA.ΔΛ / 1000;
    spac.meteo.rad_sw.e_dir .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,1) .* max(0,wd_dict["RAD_DIR"]) ./ ref_dir;
    spac.meteo.rad_sw.e_dif .= view(CACHE_CONFIG.SPECTRA.SOLAR_RAD,:,2) .* max(0,wd_dict["RAD_DIF"]) ./ ref_dif;
    spac.meteo.rad_lw = wd_dict["RAD_LW"];
    sza = solar_zenith_angle(spac.info.lat, FT(wd_dict["FDOY"]));
    spac.canopy.sun_geometry.state.sza = (wd_dict["RAD_DIR"] + wd_dict["RAD_DIF"] > 10) ? min(sza, 88.999) : sza;

    # prescribe soil water content
    swckeys = ["SWC_1", "SWC_2", "SWC_3", "SWC_4"];
    tslkeys = ["T_S_1", "T_S_2", "T_S_3", "T_S_4"];
    prescribe_soil!(spac; swcs = Tuple(min(spac.soils[i].trait.vc.Θ_SAT - 0.001, ss_dict[swckeys[i]]) for i in 1:4), t_soils = Tuple(ss_dict[tslkeys[i]] for i in 1:4));

    # prescribe leaf temperature from skin temperature
    spac.plant.memory.t_history = FT[ss_dict["T_SKN"]];
    prescribe_traits!(CACHE_CONFIG, spac; t_leaf = ss_dict["T_SKN"], t_clm = nanmean(spac.plant.memory.t_history));

    # synchronize LAI, CHL, and CI
    iday = Int(floor(wd_dict["INDEX"] / 24)) + 1;
    chl = query_griddingmachine_data(gm_dict["CHLOROPHYLL"], gm_dict["YEAR"], iday);
    ci = query_griddingmachine_data(gm_dict["CLUMPING"], gm_dict["YEAR"], iday);
    lai = query_griddingmachine_data(gm_dict["LAI"], gm_dict["YEAR"], iday);
    vcm = query_griddingmachine_data(gm_dict["VCMAX25"], gm_dict["YEAR"], iday);
    prescribe_traits!(CACHE_CONFIG, spac; ci = ci, lai = lai, vcmax = vcm, vcmax_expo = 0.3);
    prescribe_traits!(CACHE_CONFIG, spac; cab = chl, car = chl / 7);

    # initialize the spac with non-saturated soil
    initialize_states!(CACHE_CONFIG, spac);
    initialize_spac!(CACHE_CONFIG, spac);
    t_aux!(CACHE_CONFIG, spac);
    s_aux!(CACHE_CONFIG, spac);
    soil_plant_air_continuum!(CACHE_CONFIG, spac, FT(0));
    sync_state!(spac, CACHE_STATE);

    return CACHE_STATE
);

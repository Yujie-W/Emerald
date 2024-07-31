#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: add function to create data from dict
#     2024-Feb-25: add t_aux! and s_aux! functions to update the axiliary variables
#     2024-Feb-27: t_aux! and s_aux! functions moved to initialize_spac!
#     2024-Feb-28: move function to EmeraldData as it is based on GriddingMachine.jl datasets
#
#######################################################################################################################################################################################################
"""

    grid_spac(config::SPACConfiguration{FT}, gm_dict::Dict{String,Any}) where {FT}

Create a un-initialized SPAC using the data from a grid (CHL, VCMAX25, LAI, and CI are not prescribed as these changes with time), given
- `config` Configurations for SPAC
- `gm_dict` Dictionary of GriddingMachine data in a grid

"""
function grid_spac(config::SPACConfiguration{FT}, gm_dict::Dict{String,Any}) where {FT}
    #
    # TODO: add support to C4 photosynthesis
    #
    zc = max(FT(0.1), gm_dict["CANOPY_HEIGHT"]);
    spac = BulkSPAC(
                config;
                air_bounds = collect(0:21) * zc / 20,
                elevation = gm_dict["ELEVATION"],
                latitude = gm_dict["LATITUDE"],
                longitude = gm_dict["LONGITUDE"],
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                plant_zs = [-2, zc/2, zc]);
    spac.soil_bulk.trait.color = gm_dict["SOIL_COLOR"];
    @inline linear_p_soil(x) = min(1, max(eps(FT), 1 + x / 5));
    bt = BetaFunction{FT}(FUNC = linear_p_soil, PARAM_X = BetaParameterPsoil(), PARAM_Y = BetaParameterG1());
    for i in eachindex(spac.plant.leaves)
        spac.plant.leaves[i].bio.trait.lma = gm_dict["LMA"];
        #spac.plant.leaves[i].flux.trait.stomatal_model = MedlynSM{FT}(G0 = 0.005, β = bt);
        #spac.plant.leaves[i].flux.trait.stomatal_model.G1 = gm_dict["G1_MEDLYN_C3"];
        spac.plant.leaves[i].flux.trait.stomatal_model = WangSM{FT}();
    end;

    # set up SAI
    prescribe_traits!(config, spac; sai = gm_dict["SAI"]);

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

    return spac
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the spac
#     2023-Mar-13: add step to synchronize state variables into spac
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-15: make sure prescribed soil parameters are not NaN and rad is >= 0
#     2024-Feb-28: rename function to prescribe_gm_wd_data! and move it to EmeraldData
#     2024-Feb-28: make it possible to initialize spac with ss_dict
#     2024-Apr-17: update solar azimuth angle as well
#
#######################################################################################################################################################################################################
"""

    prescribe_gm_wd_data!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, ss_dict::Union{Dict{String,Any},Nothing} = nothing) where {FT}

Prescribe the SPAC with GriddingMachine and weather driver data, given
- `config` Configurations for SPAC
- `spac` SPAC to be prescribed
- `gm_dict` Dictionary of GriddingMachine data
- `wd_dict` Dictionary of weather driver data
- `ss_dict` Dictionary of initial states

"""
function prescribe_gm_wd_data!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, gm_dict::Dict{String,Any}, wd_dict::Dict{String,Any}, ss_dict::Union{Dict{String,Any},Nothing} = nothing) where {FT}
    # update environmental conditions
    for air in spac.airs
        air.state.p_air = wd_dict["P_ATM"];
        prescribe_air!(air; f_CO₂ = gm_dict["CO2"], t = wd_dict["T_AIR"], vpd = wd_dict["VPD"], wind = wd_dict["WIND"]);
    end;

    # update shortwave and longwave radiation
    ref_dir = view(config.SPECTRA.SOLAR_RAD,:,1)'  * config.SPECTRA.ΔΛ / 1000;
    ref_dif = view(config.SPECTRA.SOLAR_RAD,:,2)' * config.SPECTRA.ΔΛ / 1000;
    spac.meteo.rad_sw.e_dir .= view(config.SPECTRA.SOLAR_RAD,:,1) .* max(0,wd_dict["RAD_DIR"]) ./ ref_dir;
    spac.meteo.rad_sw.e_dif .= view(config.SPECTRA.SOLAR_RAD,:,2) .* max(0,wd_dict["RAD_DIF"]) ./ ref_dif;
    spac.meteo.rad_lw = wd_dict["RAD_LW"];
    saa = solar_azimuth_angle(spac.info.lat, FT(wd_dict["FDOY"]));
    sza = solar_zenith_angle(spac.info.lat, FT(wd_dict["FDOY"]));
    spac.canopy.sun_geometry.state.saa = saa;
    spac.canopy.sun_geometry.state.sza = (wd_dict["RAD_DIR"] + wd_dict["RAD_DIF"] > 10) ? min(sza, 88.999) : sza;

    # update t_clm to make Vcmax25 and Jmax25 TD temperature dependent
    prescribe_traits!(config, spac; t_clm = mean(spac.plant.memory.t_history));

    # synchronize LAI, CHL, and CI
    iday = Int(floor(wd_dict["INDEX"] / 24)) + 1;
    chl = query_griddingmachine_data(gm_dict["CHLOROPHYLL"], gm_dict["YEAR"], iday);
    ci = query_griddingmachine_data(gm_dict["CLUMPING"], gm_dict["YEAR"], iday);
    lai = query_griddingmachine_data(gm_dict["LAI"], gm_dict["YEAR"], iday);
    vcm = query_griddingmachine_data(gm_dict["VCMAX25"], gm_dict["YEAR"], iday);
    prescribe_traits!(config, spac; cab = chl, car = chl / 7, ci = ci, lai = lai, vcmax = vcm, vertical_expo = 0.3);

    # if ss_dict is not nothing, update soil water content and leaf temperature
    if !isnothing(ss_dict)
        # update soil water content
        swckeys = ["SWC_1", "SWC_2", "SWC_3", "SWC_4"];
        tslkeys = ["T_S_1", "T_S_2", "T_S_3", "T_S_4"];
        prescribe_soil!(spac; swcs = Tuple(min(spac.soils[i].trait.vc.Θ_SAT - 0.001, ss_dict[swckeys[i]]) for i in 1:4), t_soils = Tuple(ss_dict[tslkeys[i]] for i in 1:4));

        # prescribe leaf temperature from skin temperature
        spac.plant.memory.t_history .= FT[ss_dict["T_SKN"]];
        prescribe_traits!(config, spac; t_leaf = ss_dict["T_SKN"], t_clm = mean(spac.plant.memory.t_history));
    end;

    return nothing
end;

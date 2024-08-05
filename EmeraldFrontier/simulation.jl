#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to prescribe parameters from weather drivers
#     2023-Mar-27: prescribe T only if t_on is true, prescribe SWC only is θ_on is true
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: add a controller over rad to make sure it is >= 0
#     2023-Aug-25: make soil and leaf temperature and soil moisture optional
#     2024-Feb-22: remove state and auxil from spac struct
#     2024-Apr-17: update solar azimuth angle as well
#     2024-Jul-24: remove lai, ci, vcmax, and cab from memory (use traits instead)
# Bug fixes
#     2023-Aug-26: computed sza in the middle of a time pierod may be > 0 when cumulated radiation is > 0, set it to 88
#
#######################################################################################################################################################################################################
"""

    prescribe!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, df::NamedTuple, ind::Int; initialize_state::Bool = false) where {FT}

Prescribe traits and environmental conditions, given
- `config` `SPACConfiguration` type SPAC configuration
- `spac` `BulkSPAC` type SPAC
- `df` `NamedTuple` type weather driver
- `ind` Index of the named tuple
- `initialize_state` If true, initialize the energy state of spac

"""
function prescribe!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, df::NamedTuple, ind::Int; initialize_state::Bool = false) where {FT}
    # read the data out of dataframe row to reduce memory allocation
    df_atm::FT = df.P_ATM[ind];
    df_chl::FT = df.CHL[ind];
    df_cli::FT = df.CI[ind];
    df_co2::FT = df.CO2[ind];
    df_dif::FT = df.RAD_DIF[ind];
    df_dir::FT = df.RAD_DIR[ind];
    df_doy::FT = df.FDOY[ind];
    df_lai::FT = df.LAI[ind];
    df_lwr::FT = df.RAD_LW[ind];
    df_pcp::FT = df.PRECIP[ind];
    df_tar::FT = df.T_AIR[ind];
    df_vcm::FT = df.VCMAX25[ind];
    df_vpd::FT = df.VPD[ind];
    df_wnd::FT = df.WIND[ind];

    # prescribe the precipitation related parameters
    spac.meteo.rain = df_pcp * ρ_H₂O(FT) / M_H₂O(FT) / 3600;
    spac.meteo.t_precip = df_tar;

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    trigger_lai::Bool = !isnan(df_lai) && (df_lai != spac.canopy.structure.trait.lai);
    trigger_vcm::Bool = !isnan(df_vcm) && (df_vcm != spac.plant.leaves[end].photosystem.trait.v_cmax25);
    trigger_chl::Bool = !isnan(df_chl) && (df_chl != spac.plant.leaves[end].bio.trait.cab);
    trigger_cli::Bool = !isnan(df_cli) && (df_cli != spac.canopy.structure.trait.ci);

    if trigger_chl
        prescribe_traits!(config, spac; cab = df_chl, car = df_chl / 7);
    end;

    if trigger_vcm
        prescribe_traits!(config, spac; vcmax = df_vcm, vertical_expo = 0.3);
    end;

    if trigger_lai
        prescribe_traits!(config, spac; lai = df_lai, vertical_expo = 0.3);
    end;

    if trigger_cli
        prescribe_traits!(config, spac; ci = df_cli);
    end;

    # prescribe soil water contents and leaf temperature and initialize the spac (for first time step only)
    if initialize_state
        df_tlf::FT = df.T_LEAF[ind];
        df_ts1::FT = df.T_SOIL_1[ind];
        df_ts2::FT = df.T_SOIL_2[ind];
        df_ts3::FT = df.T_SOIL_3[ind];
        df_ts4::FT = df.T_SOIL_4[ind];
        df_sw1::FT = df.SWC_1[ind];
        df_sw2::FT = df.SWC_2[ind];
        df_sw3::FT = df.SWC_3[ind];
        df_sw4::FT = df.SWC_4[ind];

        # adjust optimum t based on the first known temperature
        @. spac.plant.memory.t_history = max(df_tar, df_tlf);
        prescribe_traits!(config, spac; t_clm = max(df_tar, df_tlf), t_leaf = max(df_tar, df_tlf));
        prescribe_soil!(spac; swcs = (df_sw1, df_sw2, df_sw3, df_sw4), t_soils = (df_ts1, df_ts2, df_ts3, df_ts4));
        initialize_spac!(config, spac);
    else
        # adjust optimum t based on 10 day moving average skin temperature
        prescribe_traits!(config, spac; t_clm = mean(spac.plant.memory.t_history));
    end;

    # update environmental conditions
    for air in spac.airs
        air.state.p_air = df_atm;
        prescribe_air!(air; f_CO₂ = df_co2, t = df_tar, vpd = df_vpd, wind = df_wnd);
    end;

    # update downward shortwave and longwave radiation
    in_dir = view(config.SPECTRA.SOLAR_RAD,:,1)' * config.SPECTRA.ΔΛ / 1000;
    in_dif = view(config.SPECTRA.SOLAR_RAD,:,2)' * config.SPECTRA.ΔΛ / 1000;
    spac.meteo.rad_sw.e_dir .= view(config.SPECTRA.SOLAR_RAD,:,1) .* max(0,df_dir) ./ in_dir;
    spac.meteo.rad_sw.e_dif .= view(config.SPECTRA.SOLAR_RAD,:,2) .* max(0,df_dif) ./ in_dif;
    spac.meteo.rad_lw = df_lwr;

    # update solar zenith angle based on the time
    saa = solar_azimuth_angle(spac.info.lat, FT(df_doy));
    sza = solar_zenith_angle(spac.info.lat, FT(df_doy));
    spac.canopy.sun_geometry.state.saa = saa;
    spac.canopy.sun_geometry.state.sza = (df_dir + df_dif > 10) ? min(sza, 88) : sza;

    # run the t_aux! and dull_aux! functions if any of the LAI, CHL, or CI changes and initialize_state is false
    if (trigger_chl || trigger_lai || trigger_cli) && !initialize_state
        t_aux!(config, spac.canopy, spac.cache);
        dull_aux!(config, spac);
    end;

    return nothing
end;


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
- `wdf` Weather driver dataframe

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

    # add the fields to store outputs
    new_df_cols = String[];
    if saving_dict["MOD_SWC"]
        for i in eachindex(spac.soils)
            push!(new_df_cols, "MOD_SWC_$i");
        end;
    end;
    if saving_dict["MOD_T_SOIL"]
        for i in eachindex(spac.soils)
            push!(new_df_cols, "MOD_T_SOIL_$i");
        end;
    end;
    if saving_dict["MOD_T_LEAF"]
        for i in eachindex(spac.plant.leaves)
            push!(new_df_cols, "MOD_T_LEAF_$i");
        end;
    end;
    if saving_dict["MOD_T_MMM"]
        push!(new_df_cols, "MOD_T_L_MAX");
        push!(new_df_cols, "MOD_T_L_MEAN");
        push!(new_df_cols, "MOD_T_L_MIN");
    end;
    if saving_dict["ΦFΦP"]
        push!(new_df_cols, "ΦF");
        push!(new_df_cols, "ΦP");
    end;
    for label in keys(saving_dict)
        if !(label in ["MOD_SWC", "MOD_T_SOIL", "MOD_T_LEAF", "MOD_T_MMM", "ΦFΦP"])
            if saving_dict[label]
                push!(new_df_cols, label);
            end;
        end;
    end;
    for label in new_df_cols
        df[!,label] .= NaN;
    end;

    # convert the DataFrame to NamedTuple
    wdf = NamedTuple{Tuple(Symbol.(names(df)))}(Tuple([df[:,n] for n in names(df)]));

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
    # read the data out of dataframe row to reduce memory allocation
    df_dif::FT = wdf.RAD_DIF[ind];
    df_dir::FT = wdf.RAD_DIR[ind];

    # prescribe parameters
    prescribe!(config, spac, wdf, ind);

    # run the model
    soil_plant_air_continuum!(config, spac, δt);
    push_t_history!(config, spac);

    # save the profiles of the soil
    if saving_dict["MOD_SWC"]
        for i in eachindex(spac.soils)
            wdf[Symbol("MOD_SWC_$i")][ind] = spac.soils[i].state.θ;
        end;
    end;
    if saving_dict["MOD_T_SOIL"]
        for i in eachindex(spac.soils)
            wdf[Symbol("MOD_T_SOIL_$i")][ind] = spac.soils[i].s_aux.t;
        end;
    end;

    # save the profiles of the leaves
    if saving_dict["MOD_T_LEAF"]
        for i in eachindex(spac.plant.leaves)
            wdf[Symbol("MOD_T_LEAF_$i")][ind] = spac.plant.leaves[i].energy.s_aux.t;
        end;
    end;
    if saving_dict["MOD_T_MMM"]
        sum_t::FT = 0;
        min_t::FT = 999;
        max_t::FT = 0;
        for l in spac.plant.leaves
            sum_t += l.energy.s_aux.t;
            min_t = min(min_t, l.energy.s_aux.t);
            max_t = max(max_t, l.energy.s_aux.t);
        end;
        wdf.MOD_T_L_MAX[ind]  = max_t;
        wdf.MOD_T_L_MEAN[ind] = sum_t / length(spac.plant.leaves);
        wdf.MOD_T_L_MIN[ind]  = min_t;
    end;

    # save the CO2 and H2O fluxes
    if saving_dict["BETA"]
        wdf.BETA[ind] = BETA(spac);
    end;
    if saving_dict["CNPP"]
        wdf.CNPP[ind] = CNPP(spac);
    end;
    if saving_dict["GPP"]
        wdf.GPP[ind] = GPP(spac);
    end;
    if saving_dict["ET_VEG"]
        wdf.ET_VEG[ind] = T_VEG(spac);
    end;

    # save the SIF (PAR and PPAR) if there is sunlight (0 otherwise)
    daytime = df_dir + df_dif >= 10;
    if saving_dict["SIF683"]
        wdf.SIF683[ind] = daytime ? TROPOMI_SIF683(config, spac) : 0;
    end;
    if saving_dict["SIF740"]
        wdf.SIF740[ind] = daytime ? TROPOMI_SIF740(config, spac) : 0;
    end;
    if saving_dict["SIF757"]
        wdf.SIF757[ind] = daytime ? OCO2_SIF759(config, spac) : 0;
    end;
    if saving_dict["SIF771"]
        wdf.SIF771[ind] = daytime ? OCO2_SIF770(config, spac) : 0;
    end;
    if saving_dict["ΣSIF"]
        wdf.ΣSIF[ind] = daytime ? ΣSIF(spac) : 0;
    end;
    if saving_dict["ΣSIF_CHL"]
        wdf.ΣSIF_CHL[ind] = daytime ? ΣSIF_CHL(config, spac) : 0;
    end;
    if saving_dict["ΣSIF_LEAF"]
        wdf.ΣSIF_LEAF[ind] = daytime ? ΣSIF_LEAF(config, spac) : 0;
    end;
    if saving_dict["PAR"]
        wdf.PAR[ind] = daytime ? PAR(config, spac) : 0;
    end;
    if saving_dict["PPAR"]
        wdf.PPAR[ind] = daytime ? PPAR(spac) : 0;
    end;

    # save the VI (and phi) if there is sunlight
    if daytime
        if saving_dict["ΦFΦP"]
            wdf.ΦF[ind],wdf.ΦP[ind] = ΦF_ΦP(spac);
        end;
        if saving_dict["NDVI"]
            wdf.NDVI[ind] = MODIS_NDVI(config, spac);
        end;
        if saving_dict["EVI"]
            wdf.EVI[ind] = MODIS_EVI(config, spac);
        end;
        if saving_dict["NIRvI"]
            wdf.NIRvI[ind] = MODIS_NIRv(config, spac);
        end;
        if saving_dict["NIRvR"]
            wdf.NIRvR[ind] = MODIS_NIRvR(config, spac);
        end;
    end;

    return nothing
);

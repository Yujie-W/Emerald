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
    trigger_cli::Bool = !isnan(df_cli) && (df_cli != spac.canopy.structure.trait.ci.ci_0);

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

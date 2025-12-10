"""

    prescribe!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, driver::NamedTuple, ind::Int; initialize_state::Bool = false) where {FT}

Prescribe traits and environmental conditions, given
- `config` `SPACConfig` type SPAC configuration
- `spac` `BulkSPAC` type SPAC
- `driver` `NamedTuple` type weather driver
- `ind` Index of the named tuple
- `initialize_state` If true, initialize the energy state of spac (forced for the first time step of a simulation)

"""
function prescribe!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, driver::NamedTuple, ind::Int; initialize_state::Bool = false) where {FT}
    # read the data out of dataframe row to reduce memory allocation
    driver_b6f::FT = driver.B6F[ind];
    driver_chl::FT = driver.CHL[ind];
    driver_cli::FT = driver.CI[ind];
    driver_co2::FT = driver.CO2[ind];
    driver_jmx::FT = driver.JMAX25[ind];
    driver_lai::FT = driver.LAI[ind];
    driver_vcm::FT = driver.VCMAX25[ind];

    driver_doy::FT = driver.FDOY[ind];
    driver_atm::FT = driver.PATM[ind];
    driver_dif::FT = driver.RAD_SW_DIF[ind];
    driver_dir::FT = driver.RAD_SW_DIR[ind];
    driver_lwr::FT = driver.RAD_LW[ind];
    driver_pcp::FT = driver.PPT[ind];
    driver_tar::FT = driver.TAIR[ind];
    driver_vpd::FT = driver.VPD[ind];
    driver_wnd::FT = driver.WIND[ind];

    # prescribe the precipitation related parameters
    # if the temperature is above T₀, it is rain, otherwise it is snow
    if driver_tar >= T₀(FT)
        spac.meteo.rain = driver_pcp * ρ_H₂O(FT) / M_H₂O(FT) / 3600;
        spac.meteo.snow = 0;
    else
        spac.meteo.rain = 0;
        spac.meteo.snow = driver_pcp * ρ_H₂O(FT) / M_H₂O(FT) / 3600;
    end;
    spac.meteo.t_precip = driver_tar;

    # vcmax, jmax, and b6f are prescribed at any time step
    prescribe_traits!(config, spac; b6f = driver_b6f, jmax = driver_jmx, vcmax = driver_vcm, vertical_expo = 0.3);

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    trigger_lai::Bool = !isnan(driver_lai) && (driver_lai != spac.canopy.structure.trait.lai);
    trigger_chl::Bool = !isnan(driver_chl) && (driver_chl != spac.plant.leaves[end].bio.trait.cab);
    trigger_cli::Bool = !isnan(driver_cli) && (driver_cli != spac.canopy.structure.trait.ci.ci_0);

    if trigger_chl
        prescribe_traits!(config, spac; cab = driver_chl, car = driver_chl / 7);
    end;

    if trigger_lai
        prescribe_traits!(config, spac; lai = driver_lai, vertical_expo = 0.3);
    end;

    if trigger_cli
        prescribe_traits!(config, spac; ci = driver_cli);
    end;

    # prescribe soil water contents and leaf temperature and initialize the spac (for first time step only)
    if initialize_state
        initialize_spac!(config, spac);
    else
        # adjust optimum t based on 10 day moving average skin temperature
        prescribe_traits!(config, spac; t_clm = mean(spac.plant.memory.t_history));
    end;

    # update environmental conditions
    for air in spac.airs
        air.state.p_air = driver_atm;
        prescribe_air!(air; f_CO₂ = driver_co2, t = driver_tar, vpd = driver_vpd, wind = driver_wnd);
    end;

    # update downward shortwave and longwave radiation
    in_dir = view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,1)' * config.CONSTANTS.SPECTRA.ΔΛ / 1000;
    in_dif = view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,2)' * config.CONSTANTS.SPECTRA.ΔΛ / 1000;
    spac.meteo.rad_sw.e_dir .= view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,1) .* max(0,driver_dir) ./ in_dir;
    spac.meteo.rad_sw.e_dif .= view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,2) .* max(0,driver_dif) ./ in_dif;
    spac.meteo.rad_lw = driver_lwr;

    # update solar zenith angle based on the time
    saa = solar_azimuth_angle(spac.info.lat, FT(driver_doy));
    sza = solar_zenith_angle(spac.info.lat, FT(driver_doy));
    spac.canopy.sun_geometry.state.saa = saa;
    spac.canopy.sun_geometry.state.sza = (driver_dir + driver_dif > 10) ? min(sza, 88) : sza;

    # run the t_aux! and dull_aux! functions if any of the LAI, CHL, or CI changes and initialize_state is false
    if (trigger_chl || trigger_lai || trigger_cli) && !initialize_state
        t_aux!(config, spac.canopy, spac.cache);
        dull_aux!(config, spac);
    end;

    return nothing
end;

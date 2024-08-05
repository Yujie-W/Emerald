#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to load all weather drivers a priori
#     2023-Mar-29: prescribe longwave radiation as well
#     2024-Feb-23: rename the function to `era5_weather_drivers` to be more specific
#     2024-Feb-29: add method to load weather drivers at a specific time index rather than loading all data
#
#######################################################################################################################################################################################################
"""

    era5_weather_drivers(wd::ERA5SingleLevelsDriver, year::Int, nx::Int)
    era5_weather_drivers(wd::ERA5SingleLevelsDriver, year::Int, nx::Int, ind::Int)

Preload weather drivers, given
- `wd` `ERA5SingleLevelsDriver` weather driver struct that stores the ERA5 info
- `year` Year of the data
- `nx` Number of grids in the 1 degree lat/lon
- `ind` Time index of the data

"""
function era5_weather_drivers end;

era5_weather_drivers(wd::ERA5SingleLevelsDriver, year::Int, nx::Int) = (
    # prescribe air layer environments and radiation
    @tinfo "Preloading weather driver for atmospheric pressure...";
    wd_p_atm = read_nc(reprocessed_file_path(wd.P_ATM[2], year, nx), wd.P_ATM[1]);
    @tinfo "Preloading weather driver for 2m air temperature...";
    wd_t_air = read_nc(reprocessed_file_path(wd.T_AIR[2], year, nx), wd.T_AIR[1]);
    @tinfo "Preloading weather driver for 2m dew temperature...";
    wd_t_dew = read_nc(reprocessed_file_path(wd.T_DEW[2], year, nx), wd.T_DEW[1]);
    @tinfo "Preloading weather driver for wind speed u...";
    wd_windu = read_nc(reprocessed_file_path(wd.WINDU[2], year, nx), wd.WINDU[1]);
    @tinfo "Preloading weather driver for wind speed v...";
    wd_windv = read_nc(reprocessed_file_path(wd.WINDV[2], year, nx), wd.WINDV[1]);
    @tinfo "Preloading weather driver for longwave radiation...";
    wd_l_all = read_nc(reprocessed_file_path(wd.L_RAD[2], year, nx), wd.L_RAD[1]);
    @tinfo "Preloading weather driver for total shortwave radiation...";
    wd_s_all = read_nc(reprocessed_file_path(wd.S_ALL[2], year, nx), wd.S_ALL[1]);
    @tinfo "Preloading weather driver for direct shortwave radiation...";
    wd_s_dir = read_nc(reprocessed_file_path(wd.S_DIR[2], year, nx), wd.S_DIR[1]);
    wd_s_dif = wd_s_all .- wd_s_dir;
    wd_vpd   = saturation_vapor_pressure.(wd_t_air) .- saturation_vapor_pressure.(wd_t_dew);
    wd_wind  = sqrt.(wd_windu .^ 2 .+ wd_windv .^ 2);

    # create a dict for the drivers
    return Dict{String,Any}(
                "YEAR"       => year,
                "RESO_SPACE" => nx,
                "P_ATM"      => wd_p_atm,
                "RAD_DIF"    => wd_s_dif,
                "RAD_DIR"    => wd_s_dir,
                "RAD_LW"     => wd_l_all,
                "T_AIR"      => wd_t_air,
                "VPD"        => wd_vpd,
                "WIND"       => wd_wind,
    );
);

era5_weather_drivers(wd::ERA5SingleLevelsDriver, year::Int, nx::Int, ind::Int) = (
    @tinfo "Load weather drivers from ERA5...";
    wd_p_atm = read_nc(reprocessed_file_path(wd.P_ATM[2], year, nx), wd.P_ATM[1], ind);
    wd_t_air = read_nc(reprocessed_file_path(wd.T_AIR[2], year, nx), wd.T_AIR[1], ind);
    wd_t_dew = read_nc(reprocessed_file_path(wd.T_DEW[2], year, nx), wd.T_DEW[1], ind);
    wd_windu = read_nc(reprocessed_file_path(wd.WINDU[2], year, nx), wd.WINDU[1], ind);
    wd_windv = read_nc(reprocessed_file_path(wd.WINDV[2], year, nx), wd.WINDV[1], ind);
    wd_l_all = read_nc(reprocessed_file_path(wd.L_RAD[2], year, nx), wd.L_RAD[1], ind);
    wd_s_all = read_nc(reprocessed_file_path(wd.S_ALL[2], year, nx), wd.S_ALL[1], ind);
    wd_s_dir = read_nc(reprocessed_file_path(wd.S_DIR[2], year, nx), wd.S_DIR[1], ind);
    wd_s_dif = wd_s_all .- wd_s_dir;
    wd_vpd   = saturation_vapor_pressure.(wd_t_air) .- saturation_vapor_pressure.(wd_t_dew);
    wd_wind  = sqrt.(wd_windu .^ 2 .+ wd_windv .^ 2);

    # create a dict for the drivers
    return Dict{String,Any}(
                "YEAR"       => year,
                "RESO_SPACE" => nx,
                "P_ATM"      => wd_p_atm,
                "RAD_DIF"    => wd_s_dif,
                "RAD_DIR"    => wd_s_dir,
                "RAD_LW"     => wd_l_all,
                "T_AIR"      => wd_t_air,
                "VPD"        => wd_vpd,
                "WIND"       => wd_wind,
    );
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-23: add function to load initial soil and skin states from ERA5
#
#######################################################################################################################################################################################################
"""

    era5_initial_states(wd::ERA5SingleLevelsDriver, year::Int, nx::Int, ind::Int)

Load initial soil and skin states from ERA5, given
- `wd` `ERA5SingleLevelsDriver` weather driver struct that stores the ERA5 info
- `year` Year of the data
- `nx` Number of grids in the 1 degree lat/lon
- `ind` Index of the data

"""
function era5_initial_states(wd::ERA5SingleLevelsDriver, year::Int, nx::Int, ind::Int)
    @tinfo "Load initial soil and skin states from ERA5...";
    states = Dict{String,Any}("YEAR" => year, "IND" => ind, "RESO_SPACE" => nx);
    states["SWC_1"] = read_nc(reprocessed_file_path(wd.SWC_1[2], year, nx), wd.SWC_1[1], ind);
    states["SWC_2"] = read_nc(reprocessed_file_path(wd.SWC_2[2], year, nx), wd.SWC_2[1], ind);
    states["SWC_3"] = read_nc(reprocessed_file_path(wd.SWC_3[2], year, nx), wd.SWC_3[1], ind);
    states["SWC_4"] = read_nc(reprocessed_file_path(wd.SWC_4[2], year, nx), wd.SWC_4[1], ind);
    states["T_S_1"] = read_nc(reprocessed_file_path(wd.T_S_1[2], year, nx), wd.T_S_1[1], ind);
    states["T_S_2"] = read_nc(reprocessed_file_path(wd.T_S_2[2], year, nx), wd.T_S_2[1], ind);
    states["T_S_3"] = read_nc(reprocessed_file_path(wd.T_S_3[2], year, nx), wd.T_S_3[1], ind);
    states["T_S_4"] = read_nc(reprocessed_file_path(wd.T_S_4[2], year, nx), wd.T_S_4[1], ind);
    states["T_SKN"] = read_nc(reprocessed_file_path(wd.T_SKN[2], year, nx), wd.T_SKN[1], ind);

    return states
end;

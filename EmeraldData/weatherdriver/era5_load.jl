#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to load all weather drivers a priori
#     2023-Mar-29: prescribe longwave radiation as well
#     2024-Feb-23: rename the function to `era5_weather_drivers` to be more specific
#
#######################################################################################################################################################################################################
"""

    era5_weather_drivers(wd::ERA5SingleLevelsDriver, year::Int, nx::Int; prescribe_soil::Bool = false)

Preload weather drivers, given
- `wd` `ERA5SingleLevelsDriver` weather driver struct that stores the ERA5 info
- `year` Year of the data
- `nx` Number of grids in the 1 degree lat/lon
- `prescribe_soil` Whether to prescribe soil water and temperature conditions, default is false (not prescribing)

"""
function era5_weather_drivers(wd::ERA5SingleLevelsDriver, year::Int, nx::Int; prescribe_soil::Bool = false)
    # prescribe air layer environments and radiation
    @tinfo "Preloading weather driver for atmospheric pressure...";
    wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(year)_$(nx)X.nc", wd.P_ATM[1]);
    @tinfo "Preloading weather driver for 2m air temperature...";
    wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(year)_$(nx)X.nc", wd.T_AIR[1]);
    @tinfo "Preloading weather driver for 2m dew temperature...";
    wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(year)_$(nx)X.nc", wd.T_DEW[1]);
    @tinfo "Preloading weather driver for wind speed u...";
    wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(year)_$(nx)X.nc", wd.WINDU[1]);
    @tinfo "Preloading weather driver for wind speed v...";
    wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(year)_$(nx)X.nc", wd.WINDV[1]);
    @tinfo "Preloading weather driver for longwave radiation...";
    wd_l_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.L_RAD[2])_SL_$(year)_$(nx)X.nc", wd.L_RAD[1]);
    @tinfo "Preloading weather driver for total shortwave radiation...";
    wd_s_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_ALL[2])_SL_$(year)_$(nx)X.nc", wd.S_ALL[1]);
    @tinfo "Preloading weather driver for direct shortwave radiation...";
    wd_s_dir = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_DIR[2])_SL_$(year)_$(nx)X.nc", wd.S_DIR[1]);
    wd_s_dif = wd_s_all .- wd_s_dir;
    wd_vpd   = saturation_vapor_pressure.(wd_t_air) .- saturation_vapor_pressure.(wd_t_dew);
    wd_wind  = sqrt.(wd_windu .^ 2 .+ wd_windv .^ 2);

    # create a dict for the drivers
    drivers = Dict{String,Any}(
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

    # prescribe soil water content
    if prescribe_soil
        @tinfo "Preloading weather drivers for soil water content at layer 1...";
        drivers["SWC_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_1[2])_SL_$(year)_$(nx)X.nc", wd.SWC_1[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 2...";
        drivers["SWC_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_2[2])_SL_$(year)_$(nx)X.nc", wd.SWC_2[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 3...";
        drivers["SWC_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_3[2])_SL_$(year)_$(nx)X.nc", wd.SWC_3[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 4...";
        drivers["SWC_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_4[2])_SL_$(year)_$(nx)X.nc", wd.SWC_4[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 1...";
        drivers["T_S_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_1[2])_SL_$(year)_$(nx)X.nc", wd.T_S_1[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 2...";
        drivers["T_S_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_2[2])_SL_$(year)_$(nx)X.nc", wd.T_S_2[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 3...";
        drivers["T_S_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_3[2])_SL_$(year)_$(nx)X.nc", wd.T_S_3[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 4...";
        drivers["T_S_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_4[2])_SL_$(year)_$(nx)X.nc", wd.T_S_4[1]);
    end;

    return drivers
end;


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
    states["SWC_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_1[2])_SL_$(year)_$(nx)X.nc", wd.SWC_1[1], ind);
    states["SWC_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_2[2])_SL_$(year)_$(nx)X.nc", wd.SWC_2[1], ind);
    states["SWC_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_3[2])_SL_$(year)_$(nx)X.nc", wd.SWC_3[1], ind);
    states["SWC_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_4[2])_SL_$(year)_$(nx)X.nc", wd.SWC_4[1], ind);
    states["T_S_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_1[2])_SL_$(year)_$(nx)X.nc", wd.T_S_1[1], ind);
    states["T_S_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_2[2])_SL_$(year)_$(nx)X.nc", wd.T_S_2[1], ind);
    states["T_S_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_3[2])_SL_$(year)_$(nx)X.nc", wd.T_S_3[1], ind);
    states["T_S_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_4[2])_SL_$(year)_$(nx)X.nc", wd.T_S_4[1], ind);
    states["T_SKN"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_SKN[2])_SL_$(year)_$(nx)X.nc", wd.T_SKN[1], ind);

    return states
end;

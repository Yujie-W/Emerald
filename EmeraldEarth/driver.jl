#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to load all weather drivers a priori
#     2023-Mar-29: prescribe longwave radiation as well
#
#######################################################################################################################################################################################################
"""

    weather_drivers(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver; leaf::Bool = true, soil::Bool = true) where {FT}

Preload weather drivers, given
- `dts` `LandDatasets` from GriddingMachine
- `wd` `ERA5SingleLevelsDriver` weather driver
- `leaf` Whether to prescribe leaf temperature from skin temperature, default is true
- `soil` Whether to prescribe soil water and temperature conditions, default is true

"""
function weather_drivers(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver; leaf::Bool = true, soil::Bool = true) where {FT}
    # prescribe air layer environments and radiation
    @tinfo "Preloading weather driver for atmospheric pressure...";
    wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1]);
    @tinfo "Preloading weather driver for 2m air temperature...";
    wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1]);
    @tinfo "Preloading weather driver for 2m dew temperature...";
    wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1]);
    @tinfo "Preloading weather driver for wind speed u...";
    wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1]);
    @tinfo "Preloading weather driver for wind speed v...";
    wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1]);
    @tinfo "Preloading weather driver for longwave radiation...";
    wd_l_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.L_RAD[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.L_RAD[1]);
    @tinfo "Preloading weather driver for total shortwave radiation...";
    wd_s_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_ALL[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_ALL[1]);
    @tinfo "Preloading weather driver for direct shortwave radiation...";
    wd_s_dir = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_DIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_DIR[1]);
    wd_s_dif = wd_s_all .- wd_s_dir;
    wd_vpd   = saturation_vapor_pressure.(wd_t_air) .- saturation_vapor_pressure.(wd_t_dew);
    wd_wind  = sqrt.(wd_windu .^ 2 .+ wd_windv .^ 2);

    # create a dict for the drivers
    drivers = Dict{String,Any}(
                "P_ATM"   => wd_p_atm,
                "RAD_DIF" => wd_s_dif,
                "RAD_DIR" => wd_s_dir,
                "RAD_LW"  => wd_l_all,
                "T_AIR"   => wd_t_air,
                "VPD"     => wd_vpd,
                "WIND"    => wd_wind,
                "YEAR"    => dts.LABELS.year,
    );

    # prescribe soil water content
    if soil
        @tinfo "Preloading weather drivers for soil water content at layer 1...";
        drivers["SWC_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_1[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 2...";
        drivers["SWC_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_2[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 3...";
        drivers["SWC_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_3[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 4...";
        drivers["SWC_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_4[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 1...";
        drivers["T_S_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_1[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 2...";
        drivers["T_S_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_2[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 3...";
        drivers["T_S_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_3[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 4...";
        drivers["T_S_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_4[1]);
    end;

    # prescribe leaf temperature from skin temperature
    if leaf
        @tinfo "Preloading weather drivers for skin temperature...";
        drivers["T_SKIN"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_SKN[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_SKN[1]);
    end;

    return drivers
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to read weather driver per grid
#     2023-Mar-13: add method to read weather driver per grid from preloaded drivers
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: add option to display process information
#
#######################################################################################################################################################################################################
"""

    wd_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; displaying::Bool = true, leaf::Bool = true, soil::Bool = true) where {FT}
    wd_grids(dts::LandDatasets{FT}, wd::Dict{String,Any}, ind::Int; displaying::Bool = true, leaf::Bool = true, soil::Bool = true) where {FT}

Prepare a matrix of weather driver data to feed SPAC, given
- `dts` `LandDatasets` from GriddingMachine
- `wd` `ERA5SingleLevelsDriver` weather driver or preloaded dictionary
- `ind` Index of data within a year
- `displaying` Whether to display information regarding process
- `leaf` Whether to prescribe leaf temperature from skin temperature, default is true
- `soil` Whether to prescribe soil water and temperature conditions, default is true

"""
function wd_grids end;

wd_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; displaying::Bool = true, leaf::Bool = true, soil::Bool = true) where {FT} = (
    if displaying
        @tinfo "Preparing a matrix of weather driver data to work on...";
    end;

    # create a matrix of GriddingMachine data
    mat_wd = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));

    # prescribe air layer environments and radiation
    wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1], ind);
    wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1], ind);
    wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1], ind);
    wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1], ind);
    wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1], ind);
    wd_l_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.L_RAD[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.L_RAD[1], ind);
    wd_s_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_ALL[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_ALL[1], ind);
    wd_s_dir = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_DIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_DIR[1], ind);
    wd_s_dif = wd_s_all .- wd_s_dir;
    wd_vpd   = saturation_vapor_pressure.(wd_t_air) .- saturation_vapor_pressure.(wd_t_dew);
    wd_wind  = sqrt.(wd_windu .^ 2 .+ wd_windv .^ 2);
    for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
        if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
            mat_wd[ilon,ilat] = Dict{String,Any}(
                        "FDOY"    => (ind - 0.5 + ((ilon - 0.5) * 360 / size(dts.t_lm,1) - 180) / 15) / 24,
                        "INDEX"   => ind,
                        "P_ATM"   => wd_p_atm[ilon,ilat],
                        "RAD_DIF" => wd_s_dif[ilon,ilat],
                        "RAD_DIR" => wd_s_dir[ilon,ilat],
                        "RAD_LW"  => wd_l_all[ilon,ilat],
                        "T_AIR"   => wd_t_air[ilon,ilat],
                        "VPD"     => wd_vpd[ilon,ilat],
                        "WIND"    => wd_wind[ilon,ilat],
                        "YEAR"    => dts.year,
            );
        end;
    end;

    # prescribe soil water content
    if soil
        wd_swc_1 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_1[1], ind);
        wd_swc_2 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_2[1], ind);
        wd_swc_3 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_3[1], ind);
        wd_swc_4 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_4[1], ind);
        wd_tsl_1 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_1[1], ind);
        wd_tsl_2 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_2[1], ind);
        wd_tsl_3 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_3[1], ind);
        wd_tsl_4 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_4[1], ind);
        for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
            if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
                mat_wd[ilon,ilat]["SWC"] = (wd_swc_1[ilon,ilat], wd_swc_2[ilon,ilat], wd_swc_3[ilon,ilat], wd_swc_4[ilon,ilat]);
                mat_wd[ilon,ilat]["T_SOIL"] = (wd_tsl_1[ilon,ilat], wd_tsl_2[ilon,ilat], wd_tsl_3[ilon,ilat], wd_tsl_4[ilon,ilat]);
            end;
        end;
    end;

    # prescribe leaf temperature from skin temperature
    if leaf
        wd_t_skn = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_SKN[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_SKN[1], ind);
        for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
            if dts.mask_spac[ilon,ilat]
                _t_skn = max(wd_t_air[ilon,ilat], wd_t_skn[ilon,ilat]);
                mat_wd[ilon,ilat]["T_SKIN"] = _t_skn;
            end;
        end;
    end;

    return mat_wd
);

wd_grids(dts::LandDatasets{FT}, wd::Dict{String,Any}, ind::Int; displaying::Bool = true, leaf::Bool = true, soil::Bool = true) where {FT} = (
    if displaying
        @tinfo "Preparing a matrix of weather driver data to work on...";
    end;

    # create a matrix of GriddingMachine data
    mat_wd = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));

    # prescribe air layer environments and radiation
    for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
        if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
            mat_wd[ilon,ilat] = Dict{String,Any}(
                        "FDOY"    => (ind - 0.5 + ((ilon - 0.5) * 360 / size(dts.t_lm,1) - 180) / 15) / 24,
                        "INDEX"   => ind,
                        "P_ATM"   => wd["P_ATM"][ilon,ilat,ind],
                        "RAD_DIF" => wd["RAD_DIF"][ilon,ilat,ind],
                        "RAD_DIR" => wd["RAD_DIR"][ilon,ilat,ind],
                        "RAD_LW"  => wd["RAD_LW"][ilon,ilat,ind],
                        "T_AIR"   => wd["T_AIR"][ilon,ilat,ind],
                        "VPD"     => wd["VPD"][ilon,ilat,ind],
                        "WIND"    => wd["WIND"][ilon,ilat,ind],
                        "YEAR"    => dts.year,
            );
        end;
    end;

    # prescribe soil water content
    if soil
        for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
            if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
                mat_wd[ilon,ilat]["SWC"] = (wd["SWC_1"][ilon,ilat,ind], wd["SWC_2"][ilon,ilat,ind], wd["SWC_3"][ilon,ilat,ind], wd["SWC_4"][ilon,ilat,ind]);
                mat_wd[ilon,ilat]["T_SOIL"] = (wd["T_S_1"][ilon,ilat,ind], wd["T_S_2"][ilon,ilat,ind], wd["T_S_3"][ilon,ilat,ind], wd["T_S_4"][ilon,ilat,ind]);
            end;
        end;
    end;

    # prescribe leaf temperature from skin temperature
    if leaf
        for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
            if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
                mat_wd[ilon,ilat]["T_SKIN"] = max(wd["T_AIR"][ilon,ilat,ind], wd["T_SKIN"][ilon,ilat,ind]);
            end;
        end;
    end;

    return mat_wd
);

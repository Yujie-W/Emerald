#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to load all weather drivers a priori
#     2023-Mar-29: prescribe longwave radiation as well
#
#######################################################################################################################################################################################################
"""

    weather_drivers(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}

Preload weather drivers, given
- `dts` `LandDatasets` from GriddingMachine
- `wd` `ERA5SingleLevelsDriver` weather driver
- `leaf` Whether to prescribe leaf temperature from skin temperature, default is true
- `soil` Whether to prescribe soil water and temperature conditions, default is true

"""
function weather_drivers(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}
    # prescribe air layer environments and radiation
    @tinfo "Preloading weather driver for atmospheric pressure...";
    _wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1]);
    @tinfo "Preloading weather driver for 2m air temperature...";
    _wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1]);
    @tinfo "Preloading weather driver for 2m dew temperature...";
    _wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1]);
    @tinfo "Preloading weather driver for wind speed u...";
    _wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1]);
    @tinfo "Preloading weather driver for wind speed v...";
    _wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1]);
    @tinfo "Preloading weather driver for longwave radiation...";
    _wd_l_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.L_RAD[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.L_RAD[1]);
    @tinfo "Preloading weather driver for total shortwave radiation...";
    _wd_s_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_ALL[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_ALL[1]);
    @tinfo "Preloading weather driver for direct shortwave radiation...";
    _wd_s_dir = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_DIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_DIR[1]);
    _wd_s_dif = _wd_s_all .- _wd_s_dir;
    _wd_vpd   = saturation_vapor_pressure.(_wd_t_air) .- saturation_vapor_pressure.(_wd_t_dew);
    _wd_wind  = sqrt.(_wd_windu .^ 2 .+ _wd_windv .^ 2);

    # create a dict for the drivers
    _drivers = Dict{String,Any}(
                "P_ATM"   => _wd_p_atm,
                "RAD_DIF" => _wd_s_dif,
                "RAD_DIR" => _wd_s_dir,
                "RAD_LW"  => _wd_l_all,
                "T_AIR"   => _wd_t_air,
                "VPD"     => _wd_vpd,
                "WIND"    => _wd_wind,
                "YEAR"    => dts.year,
    );

    # prescribe soil water content
    if soil
        @tinfo "Preloading weather drivers for soil water content at layer 1...";
        _drivers["SWC_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_1[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 2...";
        _drivers["SWC_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_2[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 3...";
        _drivers["SWC_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_3[1]);
        @tinfo "Preloading weather drivers for soil water content at layer 4...";
        _drivers["SWC_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_4[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 1...";
        _drivers["T_S_1"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_1[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 2...";
        _drivers["T_S_2"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_2[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 3...";
        _drivers["T_S_3"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_3[1]);
        @tinfo "Preloading weather drivers for soil temperature at layer 4...";
        _drivers["T_S_4"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_4[1]);
    end;

    # prescribe leaf temperature from skin temperature
    if leaf
        @tinfo "Preloading weather drivers for skin temperature...";
        _drivers["T_SKIN"] = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_SKN[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_SKN[1]);
    end;

    return _drivers
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to read weather driver per grid
#     2023-Mar-13: add method to read weather driver per grid from preloaded drivers
#     2023-Mar-29: prescribe longwave radiation as well
#
#######################################################################################################################################################################################################
"""

    wd_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}
    wd_grids(dts::LandDatasets{FT}, wd::Dict{String,Any}, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}

Prepare a matrix of weather driver data to feed SPAC, given
- `dts` `LandDatasets` from GriddingMachine
- `wd` `ERA5SingleLevelsDriver` weather driver or preloaded dictionary
- `ind` Index of data within a year
- `leaf` Whether to prescribe leaf temperature from skin temperature, default is true
- `soil` Whether to prescribe soil water and temperature conditions, default is true

"""
function wd_grids end

wd_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat} = (
    # create a matrix of GriddingMachine data
    @tinfo "Preparing a matrix of weather driver data to work on...";
    _mat_wd = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));

    # prescribe air layer environments and radiation
    _wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1], ind);
    _wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1], ind);
    _wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1], ind);
    _wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1], ind);
    _wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1], ind);
    _wd_l_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.L_RAD[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.L_RAD[1], ind);
    _wd_s_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_ALL[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_ALL[1], ind);
    _wd_s_dir = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_DIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_DIR[1], ind);
    _wd_s_dif = _wd_s_all .- _wd_s_dir;
    _wd_vpd   = saturation_vapor_pressure.(_wd_t_air) .- saturation_vapor_pressure.(_wd_t_dew);
    _wd_wind  = sqrt.(_wd_windu .^ 2 .+ _wd_windv .^ 2);
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            _mat_wd[_ilon,_ilat] = Dict{String,Any}(
                        "FDOY"    => (ind - 0.5 + ((_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180) / 15) / 24,
                        "INDEX"   => ind,
                        "P_ATM"   => _wd_p_atm[_ilon,_ilat],
                        "RAD_DIF" => _wd_s_dif[_ilon,_ilat],
                        "RAD_DIR" => _wd_s_dir[_ilon,_ilat],
                        "RAD_LW"  => _wd_l_all[_ilon,_ilat],
                        "T_AIR"   => _wd_t_air[_ilon,_ilat],
                        "VPD"     => _wd_vpd[_ilon,_ilat],
                        "WIND"    => _wd_wind[_ilon,_ilat],
                        "YEAR"    => dts.year,
            );
        end;
    end;

    # prescribe soil water content
    if soil
        _wd_swc_1 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_1[1], ind);
        _wd_swc_2 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_2[1], ind);
        _wd_swc_3 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_3[1], ind);
        _wd_swc_4 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.SWC_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.SWC_4[1], ind);
        _wd_tsl_1 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_1[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_1[1], ind);
        _wd_tsl_2 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_2[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_2[1], ind);
        _wd_tsl_3 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_3[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_3[1], ind);
        _wd_tsl_4 = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_S_4[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_S_4[1], ind);
        for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
            if dts.mask_spac[_ilon,_ilat]
                _mat_wd[_ilon,_ilat]["SWC"] = (_wd_swc_1[_ilon,_ilat], _wd_swc_2[_ilon,_ilat], _wd_swc_3[_ilon,_ilat], _wd_swc_4[_ilon,_ilat]);
                _mat_wd[_ilon,_ilat]["T_SOIL"] = (_wd_tsl_1[_ilon,_ilat], _wd_tsl_2[_ilon,_ilat], _wd_tsl_3[_ilon,_ilat], _wd_tsl_4[_ilon,_ilat]);
            end;
        end;
    end;

    # prescribe leaf temperature from skin temperature
    if leaf
        _wd_t_skn = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_SKN[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_SKN[1], ind);
        for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
            if dts.mask_spac[_ilon,_ilat]
                _t_skn = max(_wd_t_air[_ilon,_ilat], _wd_t_skn[_ilon,_ilat]);
                _mat_wd[_ilon,_ilat]["T_SKIN"] = _t_skn;
            end;
        end;
    end;

    return _mat_wd
);

wd_grids(dts::LandDatasets{FT}, wd::Dict{String,Any}, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat} = (
    # create a matrix of GriddingMachine data
    @tinfo "Preparing a matrix of weather driver data to work on...";
    _mat_wd = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));

    # prescribe air layer environments and radiation
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            _mat_wd[_ilon,_ilat] = Dict{String,Any}(
                        "FDOY"    => (ind - 0.5 + ((_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180) / 15) / 24,
                        "INDEX"   => ind,
                        "P_ATM"   => wd["P_ATM"][_ilon,_ilat,ind],
                        "RAD_DIF" => wd["RAD_DIF"][_ilon,_ilat,ind],
                        "RAD_DIR" => wd["RAD_DIR"][_ilon,_ilat,ind],
                        "RAD_LW"  => wd["RAD_LW"][_ilon,_ilat,ind],
                        "T_AIR"   => wd["T_AIR"][_ilon,_ilat,ind],
                        "VPD"     => wd["VPD"][_ilon,_ilat,ind],
                        "WIND"    => wd["WIND"][_ilon,_ilat,ind],
                        "YEAR"    => dts.year,
            );
        end;
    end;

    # prescribe soil water content
    if soil
        for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
            if dts.mask_spac[_ilon,_ilat]
                _mat_wd[_ilon,_ilat]["SWC"] = (wd["SWC_1"][_ilon,_ilat,ind], wd["SWC_2"][_ilon,_ilat,ind], wd["SWC_3"][_ilon,_ilat,ind], wd["SWC_4"][_ilon,_ilat,ind]);
                _mat_wd[_ilon,_ilat]["T_SOIL"] = (wd["T_S_1"][_ilon,_ilat,ind], wd["T_S_2"][_ilon,_ilat,ind], wd["T_S_3"][_ilon,_ilat,ind], wd["T_S_4"][_ilon,_ilat,ind]);
            end;
        end;
    end;

    # prescribe leaf temperature from skin temperature
    if leaf
        for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
            if dts.mask_spac[_ilon,_ilat]
                _mat_wd[_ilon,_ilat]["T_SKIN"] = max(wd["T_AIR"][_ilon,_ilat], wd["T_SKIN"][_ilon,_ilat]);
            end;
        end;
    end;

    return _mat_wd
);

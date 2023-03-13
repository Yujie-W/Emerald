# CLM5 settings
CLM5_PFTG = [0, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, 5.79, 5.79] .* sqrt(1000);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    spac_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}

Prepare a matrix of GriddingMachine data to feed SPAC, given
- `dts` `LandDatasets` type data struct

"""
function spac_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}
    # read some general data
    _ind_c3 = [2:14;16;17];
    _ccs = read_csv("$(@__DIR__)/../data/CO2-1Y.csv");
    _co2 = _ccs.MEAN[findfirst(_ccs.YEAR .== dts.year)];

    # create a matrix of GriddingMachine data
    @tinfo "Preparing a matrix of GriddingMachine data to work on...";
    _mat_gm = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            _pfts = dts.t_pft[_ilon,_ilat,:];
            _g = CLM5_PFTG[_ind_c3]' * _pfts[_ind_c3] / sum(_pfts[_ind_c3]);
            _g1 = isnan(_g) ? nanmean(CLM5_PFTG[_ind_c3]) : _g;
            _mat_gm[_ilon,_ilat] = Dict{String,Any}(
                        "CANOPY_HEIGHT" => dts.p_ch[_ilon,_ilat],
                        "CHLOROPHYLL"   => dts.p_chl[_ilon,_ilat,:],
                        "CLUMPING"      => dts.p_ci[_ilon,_ilat,:],
                        "CO2"           => _co2,
                        "FT"            => FT,
                        "LAI"           => dts.p_lai[_ilon,_ilat,:],
                        "LATITUDE"      => (_ilat - 0.5) * 180 / size(dts.t_lm,2) - 90,
                        "LMA"           => 1 / dts.p_sla[_ilon,_ilat] / 10,
                        "LONGITUDE"     => (_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180,
                        "MEDLYN_G1"     => _g1,
                        "SOIL_COLOR"    => min(20, max(1, Int(floor(dts.s_cc[_ilon,_ilat])))),
                        "SOIL_N"        => dts.s_n[_ilon,_ilat,:],
                        "SOIL_α"        => dts.s_α[_ilon,_ilat,:],
                        "SOIL_Θr"       => dts.s_Θr[_ilon,_ilat,:],
                        "SOIL_Θr"       => dts.s_Θr[_ilon,_ilat,:],
                        "VCMAX25"       => dts.p_vcm[_ilon,_ilat,:],
            );
        end;
    end;

    return _mat_gm
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to read weather driver per grid
#
#######################################################################################################################################################################################################
"""

    weather_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}

Prepare a matrix of weather driver data to feed SPAC, given
- `dts` `LandDatasets` from GriddingMachine
- `wd` `ERA5SingleLevelsDriver` weather driver
- `ind` Index of data within a year
- `leaf` Whether to prescribe leaf temperature from skin temperature, default is true
- `soil` Whether to prescribe soil water and temperature conditions, default is true

"""
function weather_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}
    # create a matrix of GriddingMachine data
    @tinfo "Preparing a matrix of weather driver data to work on...";
    _mat_wd = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));

    # prescribe air layer environments and radiation
    _wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1], ind);
    _wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1], ind);
    _wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1], ind);
    _wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1], ind);
    _wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1], ind);
    _wd_s_all = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_ALL[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_ALL[1], ind);
    _wd_s_dir = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.S_DIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.S_DIR[1], ind);
    _wd_s_dif = _wd_s_all .- _wd_s_dir;
    _wd_vpd   = saturation_vapor_pressure.(_wd_t_air) .- saturation_vapor_pressure.(_wd_t_dew);
    _wd_wind  = sqrt.(_wd_windu .^ 2 .+ _wd_windv .^ 2);
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            _mat_wd[_ilon,_ilat] = Dict{String,Any}(
                        "INDEX"   => ind,
                        "P_ATM"   => _wd_p_atm[_ilon,_ilat],
                        "RAD_DIF" => _wd_s_dif[_ilon,_ilat],
                        "RAD_DIR" => _wd_s_dir[_ilon,_ilat],
                        "T_AIR"   => _wd_t_air[_ilon,_ilat],
                        "UTC_DOY" => (ind .- 0.5 .+ ((_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180) / 15) ./ 24,
                        "VPD"     => _wd_vpd[_ilon,_ilat],
                        "WIND"    => _wd_wind[_ilon,_ilat],
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
end

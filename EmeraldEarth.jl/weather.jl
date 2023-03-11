#=
ERA5_LABELS = [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "mean_surface_downward_long_wave_radiation_flux",
            "mean_surface_downward_long_wave_radiation_flux_clear_sky",
            "mean_surface_direct_short_wave_radiation_flux",
            "mean_surface_direct_short_wave_radiation_flux_clear_sky",
            "mean_surface_downward_short_wave_radiation_flux",
            "mean_surface_downward_short_wave_radiation_flux_clear_sky",
            "skin_temperature",
            "total_cloud_cover",
            "total_precipitation",
ERA5_LAYERS = [
            "u10", "v10",
            "msdwlwrf", "msdwlwrfcs", "msdrswrf", "msdrswrfcs", "msdwswrf", "msdwswrfcs",
            "skt", "tcc", "tp"
ERA5_NETCDF = [
            "WIND_X", "WIND_Y",
            "RAD_LW", "RAD_LW_CS", "RAD_DIR", "RAD_DIR_CS", "RAD", "RAD_CS",
            "T_LEAF", "CLOUD", "PRECIP"
=#
#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-11: add the struct for ERA5 weather driver
#
#######################################################################################################################################################################################################
Base.@kwdef struct ERA5SingleLevelsDriver
    "Atmospheric pressure"
    P_ATM::Tuple{String,String} = ("sp", "surface_pressure")
    "Soil water content"
    SWC_1::Tuple{String,String} = ("swvl1", "volumetric_soil_water_layer_1")
    "Soil water content"
    SWC_2::Tuple{String,String} = ("swvl2", "volumetric_soil_water_layer_2")
    "Soil water content"
    SWC_3::Tuple{String,String} = ("swvl3", "volumetric_soil_water_layer_3")
    "Soil water content"
    SWC_4::Tuple{String,String} = ("swvl4", "volumetric_soil_water_layer_4")
    "Air temperature"
    T_AIR::Tuple{String,String} = ("t2m", "2m_temperature")
    "Dew temperature"
    T_DEW::Tuple{String,String} = ("d2m", "2m_dewpoint_temperature")
    "Soil temperature"
    T_S_1::Tuple{String,String} = ("stl1", "soil_temperature_level_1")
    "Soil temperature"
    T_S_2::Tuple{String,String} = ("stl2", "soil_temperature_level_2")
    "Soil temperature"
    T_S_3::Tuple{String,String} = ("stl3", "soil_temperature_level_3")
    "Soil temperature"
    T_S_4::Tuple{String,String} = ("stl4", "soil_temperature_level_4")
    "Wind speed"
    WINDU::Tuple{String,String} = ("u10", "10m_u_component_of_wind")
    "Wind speed"
    WINDV::Tuple{String,String} = ("v10", "10m_v_component_of_wind")
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to prescribe the weather drivers
#
#######################################################################################################################################################################################################
"""
"""
function prescribe! end

prescribe!(mat_spac::Matrix{Union{Nothing,MonoMLTreeSPAC}},
           dts::LandDatasets{FT},
           wd::ERA5SingleLevelsDriver,
           ind::Int;
           leaf::Bool = true,
           soil::Bool = true
) where {FT<:AbstractFloat} = (
    #@assert wd_tag in ["wd1"] "Weather driver tag not supported: $(wd_tag)";

    # prescribe air layer environments
    _wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1], ind);
    _wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1], ind);
    _wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1], ind);
    _wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1], ind);
    _wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1], ind);
    _wd_vpd   = saturation_vapor_pressure.(_wd_t_air) .- saturation_vapor_pressure.(_wd_t_dew);
    _wd_wind  = sqrt.(_wd_windu .^ 2 .+ _wd_windv .^ 2);
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        _spac = mat_spac[_ilon,_ilat];
        if _spac isa MonoMLTreeSPAC
            for _alayer in _spac.AIR
                update!(_alayer; p_CO₂ = _alayer.f_CO₂ * 1e-6 * _wd_p_atm[_ilon,_ilat], t = _wd_t_air[_ilon,_ilat], vpd = _wd_vpd[_ilon,_ilat], wind = _wd_wind[_ilon,_ilat]);
            end;
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
            _spac = mat_spac[_ilon,_ilat];
            if _spac isa MonoMLTreeSPAC
                update!(_spac;
                        swcs = (_wd_swc_1[_ilon,_ilat], _wd_swc_2[_ilon,_ilat], _wd_swc_3[_ilon,_ilat], _wd_swc_4[_ilon,_ilat]),
                        t_soils = (_wd_tsl_1[_ilon,_ilat], _wd_tsl_2[_ilon,_ilat], _wd_tsl_3[_ilon,_ilat], _wd_tsl_4[_ilon,_ilat]));
            end;
        end;
    end;

    #update!(_spac; swcs = (_df_sw1, _df_sw2, _df_sw3, _df_sw4), t_clm = _df_tmn, t_leaf = max(_df_tar, _df_tlf), t_soils = (_df_ts1, _df_ts2, _df_ts3, _df_ts4));

    return nothing
);

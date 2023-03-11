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
            "soil_temperature_level_1",
            "soil_temperature_level_2",
            "soil_temperature_level_3",
            "soil_temperature_level_4",
            "total_cloud_cover",
            "total_precipitation",
            "volumetric_soil_water_layer_1",
            "volumetric_soil_water_layer_2",
            "volumetric_soil_water_layer_3",
            "volumetric_soil_water_layer_4"];
ERA5_LAYERS = [
            "u10", "v10",
            "msdwlwrf", "msdwlwrfcs", "msdrswrf", "msdrswrfcs", "msdwswrf", "msdwswrfcs",
            "skt", "stl1", "stl2", "stl3", "stl4", "tcc", "tp", "swvl1", "swvl2", "swvl3", "swvl4"];
ERA5_NETCDF = [
            "WIND_X", "WIND_Y",
            "RAD_LW", "RAD_LW_CS", "RAD_DIR", "RAD_DIR_CS", "RAD", "RAD_CS",
            "T_LEAF", "T_SOIL_1", "T_SOIL_2", "T_SOIL_3", "T_SOIL_4", "CLOUD", "PRECIP", "SWC_1", "SWC_2", "SWC_3", "SWC_4"];
=#
Base.@kwdef struct ERA5SingleLevelsDriver
    "Atmospheric pressure"
    P_ATM::Tuple{String,String} = ("sp", "surface_pressure")
    "Air temperature"
    T_AIR::Tuple{String,String} = ("t2m", "2m_temperature")
    "Dew temperature"
    T_DEW::Tuple{String,String} = ("d2m", "2m_dewpoint_temperature")
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
           swc::Bool = true,
           t_leaf::Bool = true,
           t_soil::Bool = true
) where {FT<:AbstractFloat} = (
    #@assert wd_tag in ["wd1"] "Weather driver tag not supported: $(wd_tag)";

    # prescribe air layer environments
    @time _wd_p_atm = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.P_ATM[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.P_ATM[1], ind);
    @time _wd_t_air = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_AIR[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_AIR[1], ind);
    @time _wd_t_dew = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_DEW[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_DEW[1], ind);
    @time _wd_windu = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDU[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDU[1], ind);
    @time _wd_windv = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.WINDV[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.WINDV[1], ind);
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

    return nothing
);

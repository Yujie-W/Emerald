#=
ERA5_LABELS = [
            "mean_surface_downward_long_wave_radiation_flux",
            "mean_surface_downward_long_wave_radiation_flux_clear_sky",
            "mean_surface_direct_short_wave_radiation_flux_clear_sky",
            "mean_surface_downward_short_wave_radiation_flux_clear_sky",
            "total_cloud_cover",
            "total_precipitation",
ERA5_LAYERS = ["msdwlwrf", "msdwlwrfcs", "msdrswrfcs", "msdwswrfcs", "tcc", "tp"]
ERA5_NETCDF = ["RAD_LW", "RAD_LW_CS", "RAD_DIR_CS", "RAD_CS", "CLOUD", "PRECIP"]
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
    "Direct shortwave radiation"
    S_ALL::Tuple{String,String} = ("msdwswrf", "mean_surface_downward_short_wave_radiation_flux")
    "Direct radiation"
    S_DIR::Tuple{String,String} = ("msdrswrf", "mean_surface_direct_short_wave_radiation_flux")
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
    "Skin temperature"
    T_SKN::Tuple{String,String} = ("skt", "skin_temperature")
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
        _spac = mat_spac[_ilon,_ilat];
        if _spac isa MonoMLTreeSPAC
            for _alayer in _spac.AIR
                update!(_alayer; p_CO₂ = _alayer.f_CO₂ * 1e-6 * _wd_p_atm[_ilon,_ilat], t = _wd_t_air[_ilon,_ilat], vpd = _wd_vpd[_ilon,_ilat], wind = _wd_wind[_ilon,_ilat]);
            end;
            _in_dir = _spac.RAD_SW.e_direct' * _spac.CANOPY.WLSET.ΔΛ / 1000;
            _in_dif = _spac.RAD_SW.e_diffuse' * _spac.CANOPY.WLSET.ΔΛ / 1000;
            _spac.RAD_SW.e_direct  .= _spac.RAD_SW_REF.e_direct  .* _wd_s_dir[_ilon,_ilat] ./ _in_dir;
            _spac.RAD_SW.e_diffuse .= _spac.RAD_SW_REF.e_diffuse .* _wd_s_dif[_ilon,_ilat] ./ _in_dif;
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

    # prescribe leaf temperature from skin temperature
    if leaf
        _wd_t_skn = read_nc("$(ERA5_FOLDER)/reprocessed/$(wd.T_SKN[2])_SL_$(dts.year)_$(dts.gz)X.nc", wd.T_SKN[1], ind);
        for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
            _spac = mat_spac[_ilon,_ilat];
            if _spac isa MonoMLTreeSPAC
                _t_skn = max(_wd_t_air[_ilon,_ilat], _wd_t_skn[_ilon,_ilat]);
                push!(_spac.MEMORY.tem, _t_skn);
                if length(_spac.MEMORY.tem) > 240 deleteat!(_spac.MEMORY.tem,1) end;
                update!(_spac; t_leaf = _t_skn, t_clm = nanmean(_spac.MEMORY.tem));
            end;
        end;
    end;

    # determine whether to update LAI, CHL, and CI
    _iday = Int(floor(ind % 24)) + 1;
    _ichl = griddingmachine_data_index(size(dts.p_chl,3), dts.year, _iday);
    _icli = griddingmachine_data_index(size(dts.p_ci ,3), dts.year, _iday);
    _ilai = griddingmachine_data_index(size(dts.p_lai,3), dts.year, _iday);
    _ivcm = griddingmachine_data_index(size(dts.p_vcm,3), dts.year, _iday);
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        _spac = mat_spac[_ilon,_ilat];
        if _spac isa MonoMLTreeSPAC
            # update clumping index
            _spac.CANOPY.ci = dts.p_ci[_ilon,_ilat,_icli];
            _spac.CANOPY.Ω_A = dts.p_ci[_ilon,_ilat,_icli];

            # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
            _trigger_lai::Bool = dts.p_lai[_ilon,_ilon,_ilai] != _spac.MEMORY.lai;
            _trigger_vcm::Bool = dts.p_vcm[_ilon,_ilon,_ivcm] != _spac.MEMORY.vcm;
            _trigger_chl::Bool = dts.p_chl[_ilon,_ilon,_ichl] != _spac.MEMORY.chl;
            if _trigger_lai
                update!(_spac; lai = dts.p_lai[_ilon,_ilon,_ilai], vcmax_expo = 0.3);
                _spac.MEMORY.lai = dts.p_lai[_ilon,_ilon,_ilai];
            end;

            if _trigger_vcm
                update!(_spac; vcmax = dts.p_vcm[_ilon,_ilon,_ivcm], vcmax_expo = 0.3);
                _spac.MEMORY.vcm = dts.p_vcm[_ilon,_ilon,_ivcm];
            end;

            if _trigger_chl
                update!(_spac; cab = dts.p_chl[_ilon,_ilon,_ichl], car = dts.p_chl[_ilon,_ilon,_ichl] / 7);
                _spac.MEMORY.chl = dts.p_chl[_ilon,_ilon,_ichl];
            end;
        end;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to determine the day bounds of input GriddingMachine drivers
#
#######################################################################################################################################################################################################
"""
"""
function griddingmachine_data_index(n::Int, year::Int, d::Int)
    _bounds = [0,367];
    if n == 1
        _bounds = [0,367]
    elseif n == 12
        _bounds = isleapyear(year) ? MDAYS_LEAP : MDAYS;
    elseif n == 46
        _bounds = [collect(0:8:361); 367]
    elseif n == 52
        _bounds = [collect(0:7:361); 367]
    else
        @error "This temporal resolution is not supported: $(n)!";
    end

    return findfirst(d .<= _bounds) - 1
end

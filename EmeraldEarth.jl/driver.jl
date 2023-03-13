# ERA5 settings
ERA5_LABELS = [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "2m_dewpoint_temperature",
            "2m_temperature",
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
            "surface_pressure",
            "total_cloud_cover",
            "total_precipitation",
            "volumetric_soil_water_layer_1",
            "volumetric_soil_water_layer_2",
            "volumetric_soil_water_layer_3",
            "volumetric_soil_water_layer_4"];
ERA5_LAYERS = [
            "u10", "v10", "d2m", "t2m",
            "msdwlwrf", "msdwlwrfcs", "msdrswrf", "msdrswrfcs", "msdwswrf", "msdwswrfcs",
            "skt", "stl1", "stl2", "stl3", "stl4", "sp", "tcc", "tp", "swvl1", "swvl2", "swvl3", "swvl4"];


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-11: add the struct for ERA5 weather driver
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for ERA5 Single Levels weather driver

$(TYPEDFIELDS)

"""
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
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    fetch_ERA5_data!(year::Int; notification::Bool = false)

Fetch all the ERA5 data rquired to run this project, given
- `year` Which year of data to download
- `notification` If true, send out emails. Default is `false`

Note that ERA5 single levels data ranges from 1979 to present, and ERA5 land data ranges from 1981 to present. Be aware not to download data out of the range.

"""
function fetch_ERA5_data!(year::Int; notification::Bool = false)
    _dts = ERA5SingleLevelsHourly();
    _dir = "$(ERA5_FOLDER)/original/";
    if !isdir(_dir)
        @warn "$(_dir) does not exist, skipping...";
        return nothing
    end

    # An email will be sent out per year, comment it if you do not want
    fetch_data!(_dts, year; vars=ERA5_LABELS, folder=_dir);
    @info "Finished downloading all the datasets!";
    if notification
        send_email!("[ERA5 DATA STATUS] Downloading data for year $(year)",
                    "fluo@gps.caltech.edu",
                    "jesiner@gmail.com",
                    "ERA5 data downloading is finished for year $(year)!");
    end;

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    regrid_ERA5!(year::Int, zoom::Int = 1; notification::Bool = false)
    regrid_ERA5!(year::Int, zoom::Int, label::String, var_name::String)

Regrid the ERA5 datasets, given
- `year` Which year of data to regrid
- `zoom` The spatial resolution is `1/zoom` degree
- `notification` If true, send out emails. Default is `false`
- `label` File name label of the source NC dataset
- `var_name` Variable label in the NC dataset

To reduce memory allocation, `regrid_ERA5!` reads in the data per slice (in time index) and regrid the data per slice. What `regrid_ERA5!` does are

- determine if the output file exists
- if `true`, skip and do nothing
- if `false`
  - read the data per slice
  - replace missing values with NaNs
  - use mean of the target region to fill output matrix
  - save the regridded matrix to a new NC file

Note that ERA5 NC datasets differ from `GriddingMachine.jl` standard in that

- Latitude is from 90 to -90 in ERA5 NC dataset, and has `180N+1` elements
- Longitude is from 0 to 359.XXX in ERA5 NC dataset, and has `360N` elements

Thus, we need to regrid the dataset to ensure that the pixel orders match. For example, if the source dataset spatial resolution is `4X` and the target spatial resolution is `1X`, we need take the
    average of data in the pixel of 120E-121E (the steps are 120, 120.25, 120.5, 120.75, 121), 30N-31N (the steps are 30, 30.25, 30.5, 30.75, 31), namely 25 elements in total. The special case is
    when the longitude ranges from 179E to 180E, and in this case, we need to include the -180E slice. Otherwise, the mean longitude is 179.375 rather than 179.5. The general function is

"""
function regrid_ERA5! end

regrid_ERA5!(year::Int, zoom::Int = 1; notification::Bool = false) = (
    regrid_ERA5!.(year, zoom, ERA5_LABELS, ERA5_LAYERS);
    @info "Finished regridding all the datasets!";
    if notification
        send_email!("[ERA5 DATA STATUS] Regridding data for year $(year)",
                    "fluo@gps.caltech.edu",
                    "jesiner@gmail.com",
                    "ERA5 data regridding is finished for year $(year)!");
    end;

    return nothing;
);

regrid_ERA5!(year::Int, zoom::Int, label::String, var_name::String) = (
    _file_in  = "$(ERA5_FOLDER)/original/$(label)_SL_$(year).nc";
    _file_out = "$(ERA5_FOLDER)/reprocessed/$(label)_SL_$(year)_$(zoom)X.nc";

    # if file exists already, skip
    if isfile(_file_out)
        @info "File $(_file_out) already exists!";
        return nothing;
    end;

    # read the file per slice
    @info "Reading and regridding file $(_file_in) per time slice...";
    _time = read_nc(_file_in, "time");
    _matn = zeros((360*zoom,180*zoom,length(_time))) .* NaN;
    @showprogress for _itim in eachindex(_time)
        _mati = read_nc(_file_in, var_name, _itim);
        _z_ss = Int(size(_mati,1) / 360);
        _dxy  = Int(_z_ss/zoom);
        _nvar = replace(_mati, missing=>NaN);
        for _ilon in axes(_matn,1), _ilat in axes(_matn,2)
            # from +0 to +180
            if size(_matn,1)/2 < _ilon <= size(_matn,1)
                __ilon = _ilon - Int(size(_matn,1)/2);
                _sub   = _nvar[_dxy*(__ilon-1)+1:_dxy*__ilon+1, _dxy*(_ilat-1)+1:_dxy*_ilat+1];
            # from -180 to -0.NNN
            elseif _ilon < size(_matn,1)/2
                __ilon = _ilon + Int(size(_matn,1)/2);
                _sub   = _nvar[_dxy*(__ilon-1)+1:_dxy*__ilon+1, _dxy*(_ilat-1)+1:_dxy*_ilat+1];
            # -0
            else
                __ilon = size(_matn,1);
                _sub   = [collect(_nvar[_dxy*(__ilon-1)+1:_dxy*__ilon, _dxy*(_ilat-1)+1:_dxy*_ilat+1]); collect(_nvar[1:1, _dxy*(_ilat-1)+1:_dxy*_ilat+1])];
            end;
            # for soil water content, replace with NaN if the SWC <= 0.01
            if var_name in ["SWC_1", "SWC_2", "SWC_3", "SWC_4"]
                _sub[_sub .<= 0] .= NaN;
            end;
            _matn[_ilon,size(_matn,2)+1-_ilat,_itim] = nanmean(_sub);
        end;
    end;

    # save the regridded dataset
    @info "Saving regridded dataset to $(_file_out)...";
    _attr = Dict(var_name => label * "_SL", "unit" => "Same as $(_file_in)");
    save_nc!(_file_out, var_name, _matn, _attr);

    return nothing;
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to read weather driver per grid
#
#######################################################################################################################################################################################################
"""

    wd_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}

Prepare a matrix of weather driver data to feed SPAC, given
- `dts` `LandDatasets` from GriddingMachine
- `wd` `ERA5SingleLevelsDriver` weather driver
- `ind` Index of data within a year
- `leaf` Whether to prescribe leaf temperature from skin temperature, default is true
- `soil` Whether to prescribe soil water and temperature conditions, default is true

"""
function wd_grids(dts::LandDatasets{FT}, wd::ERA5SingleLevelsDriver, ind::Int; leaf::Bool = true, soil::Bool = true) where {FT<:AbstractFloat}
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
                        "FDOY"    => (ind - 0.5 + ((_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180) / 15) / 24,
                        "INDEX"   => ind,
                        "P_ATM"   => _wd_p_atm[_ilon,_ilat],
                        "RAD_DIF" => _wd_s_dif[_ilon,_ilat],
                        "RAD_DIR" => _wd_s_dir[_ilon,_ilat],
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
end

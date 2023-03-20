module ERA5

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Fetcher: ERA5SingleLevelsHourly, fetch_data!

using ..EmeraldCore.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldIO.Netcdf: append_nc!, read_nc, save_nc!, varname_nc
using ..EmeraldMath.Stats: nanmean
using ..EmeraldUtility.Email: send_email!


# CliMA Land settings
DRIVER_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/drivers";

# ERA5 settings
ERA5_FOLDER = "/home/wyujie/DATASERVER/reanalysis/ERA5/SingleLevels";
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
ERA5_NETCDF = [
            "WIND_X", "WIND_Y", "T_DEW", "T_AIR",
            "RAD_LW", "RAD_LW_CS", "RAD_DIR", "RAD_DIR_CS", "RAD", "RAD_CS",
            "T_LEAF", "T_SOIL_1", "T_SOIL_2", "T_SOIL_3", "T_SOIL_4", "P_ATM", "CLOUD", "PRECIP", "SWC_1", "SWC_2", "SWC_3", "SWC_4"];


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
#     2023-Mar-20: move function from ClimaLand-0.2
#     2023-Mar-20: add code to read weather driver for a single grid
#
#######################################################################################################################################################################################################
"""

    weather_driver_file(wd_tag::String, dict::Dict{String,Any}; displaying::Bool = true, appending::Bool = false)

Return the input weather driver file path, given
- `wd_tag` Weather driver version tag
- `dict` Dictionary that store grid information
- `displaying` If true, display information about the NetCDF file
- `appending` If true, always check whether there are new fields to add

"""
function weather_driver_file(wd_tag::String, dict::Dict{String,Any}; displaying::Bool = true, appending::Bool = false)
    # which index of data to read
    _gz      = dict["RESO_SPACE"]
    _lat_ind = dict["LAT_INDEX"];
    _lon     = dict["LONGITUDE"];
    _lon_ind = dict["LON_INDEX"];
    _year    = dict["YEAR"];

    # folders that stores the input data
    @assert isdir(DRIVER_FOLDER) "Weather driver folder $(DRIVER_FOLDER) does not exist...";
    _nc_name = "weather_driver_$(wd_tag)_$(_year)_$(_lat_ind)_$(_lon_ind)_$(_gz).nc";
    _nc_path = "$(DRIVER_FOLDER)/$(_year)/$(_nc_name)";

    # if file exists and appending is false
    if isfile(_nc_path) && !appending
        if displaying
            @info "$(_nc_path) exists, doing nothing...";
        end;

        return _nc_path
    end;

    # if file exists and appending is true
    if isfile(_nc_path) && !appending
        if displaying
            @info "$(_nc_path) exists, adding new fields...";
        end;

        _existed_varnames = varname_nc(_nc_path);
        for _i in eachindex(ERA5_LABELS)
            if !(ERA5_NETCDF[_i] in _existed_varnames)
                _nc_var = "$(ERA5_FOLDER)/reprocessed/$(ERA5_LABELS[_i])_SL_$(_year)_$(_gz).nc";
                _nc_vec = read_nc(_nc_var, ERA5_LAYERS[_i], _lon_ind, _lat_ind);
                append_nc!(_nc_path, ERA5_NETCDF[_i], _nc_vec, Dict{String,String}("longname" => ERA5_LABELS[_i] * "_SL"), ["ind"]);
            end;
        end;

        return _nc_path
    end;

    # if file does not exist
    if displaying
        @info "$(_nc_path) does not exist, generating file now...";
    end;

    # function to create DataFrame columns
    @inline function add_col!(df::DataFrame, label::String, layer::String, var_name::String)
        _nc_var = "$(ERA5_FOLDER)/reprocessed/$(label)_SL_$(_year)_$(_gz).nc";
        df[!, var_name] = read_nc(_nc_var, layer, _lon_ind, _lat_ind);

        return _nc_path;
    end;

    # add data into DataFrame
    _tz = _lon / 15;
    _df = DataFrame();
    add_col!.([_df], ERA5_LABELS, ERA5_LAYERS, ERA5_NETCDF);

    _df[!,"FDOY"   ] = (collect(eachindex(_df.P_ATM)) .- 0.5 .+ _tz) ./ 24;
    _df[!,"WIND"   ] = sqrt.( _df.WIND_X .^ 2 .+ _df.WIND_Y .^2 );
    _df[!,"RAD_DIF"] = _df.RAD .- _df.RAD_DIR;
    _df[!,"VPD"    ] = saturation_vapor_pressure.(_df.T_AIR) .- saturation_vapor_pressure.(_df.T_DEW);
    _var_attrs = Dict{String,String}[[Dict{String,String}("longname" => ERA5_LABELS[_i] * "_SL") for _i in eachindex(ERA5_LABELS)];
                                     Dict{String,String}("longname" => "Day of year");
                                     Dict{String,String}("longname" => "Wind speed");
                                     Dict{String,String}("longname" => "Diffuse radiation");
                                     Dict{String,String}("longname" => "Vapor pressure deficit")];
    save_nc!(_nc_path, _df, [ERA5_NETCDF; "FDOY"; "WIND"; "RAD_DIF"; "VPD"], _var_attrs);

    return _nc_path
end


end # module

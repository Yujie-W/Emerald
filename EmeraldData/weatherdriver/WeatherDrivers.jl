module WeatherDrivers

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Fetcher: fetch_data!
using NetcdfIO: append_nc!, read_nc, save_nc!, varname_nc

using ..EmeraldLand.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldMath.Stats: nanmean
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Email: send_email!


include("era5_regrid.jl");
include("era5_setting.jl");
include("era5_type.jl");

include("era5_load.jl");

include("parser.jl");


# CliMA Land settings
DRIVER_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/drivers";


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: move function from ClimaLand-0.2
#     2023-Mar-20: add code to read weather driver for a single grid
#     2023-May-11: fix a typo in appending
#
#######################################################################################################################################################################################################
"""

    weather_driver_file(wd_tag::String, dict::Dict{String,Any}; appending::Bool = false)

Return the input weather driver file path and name, given
- `wd_tag` Weather driver version tag
- `dict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add

"""
function weather_driver_file(wd_tag::String, dict::Dict{String,Any}; appending::Bool = false)
    # which index of data to read
    nx       = dict["RESO_SPACE"]
    _lat_ind = dict["LAT_INDEX"];
    _lon     = dict["LONGITUDE"];
    _lon_ind = dict["LON_INDEX"];
    _year    = dict["YEAR"];
    msg_lvl  = dict["MESSAGE_LEVEL"];

    # folders that stores the input data
    @assert isdir(DRIVER_FOLDER) "Weather driver folder $(DRIVER_FOLDER) does not exist...";
    _nc_name = "weather_driver_$(wd_tag)_$(_year)_$(_lat_ind)_$(_lon_ind)_$(nx)X.nc";
    _nc_path = "$(DRIVER_FOLDER)/$(_year)/$(_nc_name)";

    # if file exists and appending is false
    if isfile(_nc_path) && !appending
        if msg_lvl == 2
            @info "$(_nc_path) exists, doing nothing...";
        end;

        return _nc_path, _nc_name
    end;

    # if file exists and appending is true
    if isfile(_nc_path) && appending
        if msg_lvl == 2
            @info "$(_nc_path) exists, adding new fields...";
        end;

        _existed_varnames = varname_nc(_nc_path);
        for i in eachindex(ERA5_LABELS)
            if !(ERA5_NETCDF[i] in _existed_varnames)
                _nc_var = "$(ERA5_FOLDER)/reprocessed/$(ERA5_LABELS[i])_SL_$(_year)_$(nx)X.nc";
                _nc_vec = read_nc(_nc_var, ERA5_LAYERS[i], _lon_ind, _lat_ind);
                append_nc!(_nc_path, ERA5_NETCDF[i], _nc_vec, Dict{String,String}("longname" => ERA5_LABELS[i] * "_SL"), ["ind"]);
            end;
        end;

        return _nc_path, _nc_name
    end;

    # if file does not exist
    if msg_lvl == 2
        @info "$(_nc_path) does not exist, generating file now...";
    end;

    # function to create DataFrame columns
    @inline function add_col!(df::DataFrame, label::String, layer::String, var_name::String)
        _nc_var = "$(ERA5_FOLDER)/reprocessed/$(label)_SL_$(_year)_$(nx)X.nc";
        df[!, var_name] = read_nc(_nc_var, layer, _lon_ind, _lat_ind);

        return _nc_path, _nc_name
    end;

    # add data into DataFrame
    _tz = _lon / 15;
    _df = DataFrame();
    add_col!.([_df], ERA5_LABELS, ERA5_LAYERS, ERA5_NETCDF);

    _df[!,"FDOY"   ] = (collect(eachindex(_df.P_ATM)) .- 0.5 .+ _tz) ./ 24;
    _df[!,"WIND"   ] = sqrt.( _df.WIND_X .^ 2 .+ _df.WIND_Y .^2 );
    _df[!,"RAD_DIF"] = _df.RAD .- _df.RAD_DIR;
    _df[!,"VPD"    ] = saturation_vapor_pressure.(_df.T_AIR) .- saturation_vapor_pressure.(_df.T_DEW);
    _var_attrs = Dict{String,String}[[Dict{String,String}("longname" => ERA5_LABELS[i] * "_SL") for i in eachindex(ERA5_LABELS)];
                                      Dict{String,String}("longname" => "Day of year");
                                      Dict{String,String}("longname" => "Wind speed");
                                      Dict{String,String}("longname" => "Diffuse radiation");
                                      Dict{String,String}("longname" => "Vapor pressure deficit")];
    save_nc!(_nc_path, _df, [ERA5_NETCDF; "FDOY"; "WIND"; "RAD_DIF"; "VPD"], _var_attrs);

    return _nc_path, _nc_name
end;


end; # module

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
    lat_ind  = dict["LAT_INDEX"];
    lon      = dict["LONGITUDE"];
    lon_ind  = dict["LON_INDEX"];
    msg_lvl  = dict["MESSAGE_LEVEL"];
    nx       = dict["RESO_SPACE"]
    year     = dict["YEAR"];

    # folders that stores the input data
    @assert isdir(DRIVER_FOLDER) "Weather driver folder $(DRIVER_FOLDER) does not exist...";
    nc_name = "weather_driver_$(wd_tag)_$(year)_$(lat_ind)_$(lon_ind)_$(nx)X.nc";
    nc_path = "$(DRIVER_FOLDER)/$(year)/$(nc_name)";

    # if file exists and appending is false
    if isfile(nc_path) && !appending
        if msg_lvl == 2
            @info "$(nc_path) exists, doing nothing...";
        end;

        return nc_path, nc_name
    end;

    # if file exists and appending is true
    if isfile(nc_path) && appending
        if msg_lvl == 2
            @info "$(nc_path) exists, adding new fields...";
        end;

        existed_varnames = varname_nc(nc_path);
        for i in eachindex(ERA5_LABELS)
            if !(ERA5_NETCDF[i] in existed_varnames)
                nc_var = "$(ERA5_FOLDER)/reprocessed/$(ERA5_LABELS[i])_$(year)_$(nx)X.nc";
                nc_vec = read_nc(nc_var, ERA5_LAYERS[i], lon_ind, lat_ind);
                append_nc!(nc_path, ERA5_NETCDF[i], nc_vec, Dict{String,String}("longname" => ERA5_LABELS[i]), ["ind"]);
            end;
        end;

        return nc_path, nc_name
    end;

    # if file does not exist
    if msg_lvl == 2
        @info "$(nc_path) does not exist, generating file now...";
    end;

    # function to create DataFrame columns
    @inline function add_col!(df::DataFrame, label::String, layer::String, var_name::String)
        nc_var = "$(ERA5_FOLDER)/reprocessed/$(label)_$(year)_$(nx)X.nc";
        df[!, var_name] = read_nc(nc_var, layer, lon_ind, lat_ind);

        return nc_path, nc_name
    end;

    # add data into DataFrame
    tz = lon / 15;
    df = DataFrame();
    add_col!.([df], ERA5_LABELS, ERA5_LAYERS, ERA5_NETCDF);

    df[!,"FDOY"   ] = (collect(eachindex(df.P_ATM)) .- 0.5 .+ tz) ./ 24;
    df[!,"WIND"   ] = sqrt.( df.WIND_X .^ 2 .+ df.WIND_Y .^2 );
    df[!,"RAD_DIF"] = df.RAD .- df.RAD_DIR;
    df[!,"VPD"    ] = saturation_vapor_pressure.(df.T_AIR) .- saturation_vapor_pressure.(df.T_DEW);
    var_attrs = Dict{String,String}[[Dict{String,String}("longname" => ERA5_LABELS[i]) for i in eachindex(ERA5_LABELS)];
                                     Dict{String,String}("longname" => "Day of year");
                                     Dict{String,String}("longname" => "Wind speed");
                                     Dict{String,String}("longname" => "Diffuse radiation");
                                     Dict{String,String}("longname" => "Vapor pressure deficit")];
    save_nc!(nc_path, df, [ERA5_NETCDF; "FDOY"; "WIND"; "RAD_DIF"; "VPD"], var_attrs);

    return nc_path, nc_name
end;

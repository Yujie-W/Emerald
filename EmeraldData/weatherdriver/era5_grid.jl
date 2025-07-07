#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: move function from ClimaLand-0.2
#     2023-Mar-20: add code to read weather driver for a single grid
#     2023-May-11: fix a typo in appending
#     2024-Mar-07: rename function to era5_weather_driver_file using the ERA5SingleLevelsDriver struct
#
#######################################################################################################################################################################################################
"""

    era5_weather_driver_file(wd::ERA5SingleLevelsDriver, gm_dict::Dict{String,Any}; appending::Bool = false)

Return the ERA5 weather driver file path and name, given
- `wd` `ERA5SingleLevelsDriver` weather driver struct that stores the ERA5 info
- `dict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add

"""
function era5_weather_driver_file(wd::ERA5SingleLevelsDriver, gm_dict::Dict{String,Any}; appending::Bool = false)
    # which index of data to read
    lat_ind  = gm_dict["LAT_INDEX"];
    lon      = gm_dict["LONGITUDE"];
    lon_ind  = gm_dict["LON_INDEX"];
    msg_lvl  = gm_dict["MESSAGE_LEVEL"];

    # folders that stores the input data
    @assert isdir(LAND_DRIVER) "Weather driver folder $(LAND_DRIVER) does not exist...";
    nc_path = grid_file_path(gm_dict);

    # if file exists and appending is false
    if isfile(nc_path) && !appending
        if msg_lvl == 2
            @info "$(nc_path) exists, doing nothing...";
        end;

        return nc_path
    end;

    # if file exists and appending is true
    if isfile(nc_path) && appending
        if msg_lvl == 2
            @info "$(nc_path) exists, adding new fields...";
        end;

        existed_varnames = varname_nc(nc_path);
        for fn in fieldnames(ERA5SingleLevelsDriver)
            varfn = getfield(wd, fn);
            if !(varfn[3] in existed_varnames)
                nc_var = reprocessed_file_path(gm_dict, varfn[2]);
                nc_vec = read_nc(nc_var, varfn[1], lon_ind, lat_ind);
                append_nc!(nc_path, varfn[3], nc_vec, Dict{String,String}("longname" => varfn[2]), ["ind"]);
            end;
        end;

        return nc_path
    end;

    # if file does not exist
    if msg_lvl == 2
        @info "$(nc_path) does not exist, generating file now...";
    end;

    # function to create DataFrame columns
    @inline function add_col!(df::DataFrame, label::String, layer::String, var_name::String)
        nc_var = reprocessed_file_path(gm_dict, label);
        df[!, var_name] = read_nc(nc_var, layer, lon_ind, lat_ind);

        return nc_path
    end;

    # add data into DataFrame
    tz = lon / 15;
    df = DataFrame();
    for fn in fieldnames(ERA5SingleLevelsDriver)
        varfn = getfield(wd, fn);
        add_col!(df, varfn[2], varfn[1], varfn[3]);
    end;

    df[!,"FDOY"   ] = (collect(eachindex(df.P_ATM)) .- 0.5 .+ tz) ./ 24;
    df[!,"WIND"   ] = sqrt.( df.WIND_X .^ 2 .+ df.WIND_Y .^2 );
    df[!,"RAD_DIF"] = df.RAD .- df.RAD_DIR;
    df[!,"VPD"    ] = saturation_vapor_pressure.(df.T_AIR) .- saturation_vapor_pressure.(df.T_DEW);
    var_labels = [getfield(wd, fn)[2] for fn in fieldnames(ERA5SingleLevelsDriver)];
    var_dflabs = [getfield(wd, fn)[3] for fn in fieldnames(ERA5SingleLevelsDriver)];
    var_attrs = Dict{String,String}[[Dict{String,String}("longname" => label) for label in var_labels];
                                     Dict{String,String}("longname" => "Day of year");
                                     Dict{String,String}("longname" => "Wind speed");
                                     Dict{String,String}("longname" => "Diffuse radiation");
                                     Dict{String,String}("longname" => "Vapor pressure deficit")];
    save_nc!(nc_path, df, [var_dflabs; "FDOY"; "WIND"; "RAD_DIF"; "VPD"], var_attrs);

    return nc_path
end;

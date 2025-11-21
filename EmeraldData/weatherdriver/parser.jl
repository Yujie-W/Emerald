#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-23: add function to preload weather drivers based on the tag
#     2024-Feb-29: add method to load weather drivers at a specific time index rather than loading all data
#
#######################################################################################################################################################################################################
"""

    preloaded_weather_drivers(wd_tag::String, year::Int, nx::Int)
    preloaded_weather_drivers(wd_tag::String, year::Int, nx::Int, ind::Int)

Preload weather drivers, given
- `wd_tag` Weather driver version tag
- `year` Year of the data
- `nx` Number of grids in the 1 degree lat/lon
- `ind` Time index of the data

"""
function preloaded_weather_drivers end;

preloaded_weather_drivers(wd_tag::String, year::Int, nx::Int) = (
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    if wd_tag == "wd1"
        return era5_weather_drivers(ERA5SingleLevelsDriver(), year, nx)
    end;
);

preloaded_weather_drivers(wd_tag::String, year::Int, nx::Int, ind::Int) = (
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    if wd_tag == "wd1"
        return era5_weather_drivers(ERA5SingleLevelsDriver(), year, nx, ind)
    end;
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: move function from ClimaLand-0.2
#     2023-Mar-28: add simulated swc and temperatures into dataframe so as to output
#     2023-Aug-25: move method to interpolate data to EmeraldMath.jl
#     2024-Mar-07: move function from EmeraldFrontier to EmeraldData (to use with EmeraldFrontier and ClimaLandPRO)
#     2025-Sep-01: add new methods to read from an existing weather driver file
#
#######################################################################################################################################################################################################
"""

    grid_weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false)

Prepare weather driver dataframe in a grid to feed SPAC, given
- `wd_tag` Weather driver version tag
- `gm_dict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add

"""
function grid_weather_driver end;

grid_weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false) = (
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    if wd_tag == "wd1"
        nc_wd = era5_weather_driver_file(ERA5SingleLevelsDriver(), gm_dict; appending = appending);

        return grid_weather_driver(wd_tag, gm_dict, nc_wd)
    end;

    return error("Weather driver tag $(wd_tag) is not supported...")
);

grid_weather_driver(wd_tag::String, gm_dict::Dict{String,Any}, nc_path::String) = (
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    df_wd = read_nc(nc_path);

    # interpolate the data to a new resolution
    df_wd[!,"CO2"    ] .= resample(gm_dict["CO2"        ], "1H", gm_dict["YEAR"]);
    df_wd[!,"CHL"    ] .= resample(gm_dict["CHLOROPHYLL"], "1H", gm_dict["YEAR"]);
    df_wd[!,"CI"     ] .= resample(gm_dict["CLUMPING"   ], "1H", gm_dict["YEAR"]);
    df_wd[!,"LAI"    ] .= resample(gm_dict["LAI"        ], "1H", gm_dict["YEAR"]);
    df_wd[!,"VCMAX25"] .= resample(gm_dict["VCMAX25"    ], "1H", gm_dict["YEAR"]);

    return df_wd
);

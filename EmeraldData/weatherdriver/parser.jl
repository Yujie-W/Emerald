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
#     2024-Feb-23: add function to initial soil and skin states based on the tag
#
#######################################################################################################################################################################################################
"""

    initial_soil_skin_states(wd_tag::String, year::Int, nx::Int, ind::Int)

Initial soil and skin states, given
- `wd_tag` Weather driver version tag
- `year` Year of the data
- `nx` Number of grids in the 1 degree lat/lon
- `ind` Index of the data

"""
function initial_soil_skin_states(wd_tag::String, year::Int, nx::Int, ind::Int)
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    if wd_tag == "wd1"
        return era5_initial_states(ERA5SingleLevelsDriver(), year, nx, ind)
    end;
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: move function from ClimaLand-0.2
#     2023-Mar-28: add simulated swc and temperatures into dataframe so as to output
#     2023-Aug-25: move method to interpolate data to EmeraldMath.jl
#     2024-Mar-07: move function from EmeraldFrontier to EmeraldData (to use with EmeraldFrontier and ClimaLandPRO)
#
#######################################################################################################################################################################################################
"""

    grid_weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false)

Prepare weather driver dataframe in a grid to feed SPAC, given
- `wd_tag` Weather driver version tag
- `gm_dict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add

"""
function grid_weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false)
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    if wd_tag == "wd1"
        nc_wd = era5_weather_driver_file(ERA5SingleLevelsDriver(), gm_dict; appending = appending);
        df_wd = read_nc(nc_wd);

        # interpolate the data to a new resolution
        df_wd[!,"CO2"    ] .= interpolate_data(gm_dict["CO2"], gm_dict["YEAR"]; out_reso = "1H");
        df_wd[!,"CHL"    ] .= interpolate_data(gm_dict["CHLOROPHYLL"], gm_dict["YEAR"]; out_reso = "1H");
        df_wd[!,"CI"     ] .= interpolate_data(gm_dict["CLUMPING"], gm_dict["YEAR"]; out_reso = "1H");
        df_wd[!,"LAI"    ] .= interpolate_data(gm_dict["LAI"], gm_dict["YEAR"]; out_reso = "1H");
        df_wd[!,"VCMAX25"] .= interpolate_data(gm_dict["VCMAX25"], gm_dict["YEAR"]; out_reso = "1H");
    end;

    return df_wd
end;

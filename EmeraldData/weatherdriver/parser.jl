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

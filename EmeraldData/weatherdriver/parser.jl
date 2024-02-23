#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-23: add function to preload weather drivers based on the tag
#
#######################################################################################################################################################################################################
"""

    preloaded_weather_drivers(wd_tag::String, year::Int, nx::Int; prescribe_soil::Bool = false)

Preload weather drivers, given
- `wd_tag` Weather driver version tag
- `year` Year of the data
- `nx` Number of grids in the 1 degree lat/lon
- `prescribe_soil` Whether to prescribe soil water and temperature conditions, default is false (not prescribing)

"""
function preloaded_weather_drivers(wd_tag::String, year::Int, nx::Int; prescribe_soil::Bool = false)
    @assert wd_tag in ["wd1"] "Weather driver tag $(wd_tag) is not supported...";

    # wd1 is the ERA5 single level driver
    if wd_tag == "wd1"
        return era5_weather_drivers(ERA5SingleLevelsDriver(), year, nx; prescribe_soil = prescribe_soil)
    end;
end;


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

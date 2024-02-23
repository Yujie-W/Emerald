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
#     2024-Feb-23: add function to extract a slice of weather data from preloaded drivers
#
#######################################################################################################################################################################################################
"""

    weather_drivers_snapshot(wds::Dict{String,Any}, ind::Int)

Extract a slice of weather data from preloaded drivers, given
- `wds` Preloaded weather drivers
- `ind` Index of the data

"""
function weather_drivers_snapshot(wds::Dict{String,Any}, ind::Int)
    @tinfo "Extract a slice of weather data from preloaded drivers...";

    # create a matrix of GriddingMachine data
    wd_keys = [k for k in keys(wds) if k != "RESO_SPACE" && k != "YEAR"];
    nx = wds["RESO_SPACE"];
    mat_wd = Matrix{Dict{String,Any}}(undef, 360nx, 180nx);
    for ilon in axes(mat_wd,1), ilat in axes(mat_wd,2)
        dict = Dict{String,Any}("YEAR" => wds["YEAR"], "RESO_SPACE" => wds["RESO_SPACE"]);
        for k in wd_keys
            dict[k] = wds[k][ilon,ilat,ind];
        end;
        mat_wd[ilon,ilat] = dict;
    end;

    return mat_wd
end;

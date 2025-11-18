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
- `ind` Index of the hourly data

"""
function weather_drivers_snapshot(wds::Dict{String,Any}, ind::Int)
    pretty_display!("Extracting a slice matrix of weather data from preloaded drivers...", "tinfo");

    # create a matrix of weather drivers
    wd_keys = [k for k in keys(wds) if !(k in ["RESO_SPACE", "YEAR", "IND"])];
    nx = wds["RESO_SPACE"];
    reso = 1 / nx;
    mat_wd = Matrix{Dict{String,Any}}(undef, 360nx, 180nx);
    for ilon in axes(mat_wd,1), ilat in axes(mat_wd,2)
        dict = Dict{String,Any}(
                    "YEAR" => wds["YEAR"],
                    "RESO_SPACE" => wds["RESO_SPACE"],
                    "INDEX" => ind,
                    "FDOY" => ((-180 + (ilon - 0.5) * reso) / 15 + ind - 0.5) / 24);
        for k in wd_keys
            dict[k] = wds[k][ilon,ilat,ind];
        end;
        mat_wd[ilon,ilat] = dict;
    end;

    return mat_wd
end;

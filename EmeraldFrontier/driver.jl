#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: move function from ClimaLand-0.2
#     2023-Mar-28: add simulated swc and temperatures into dataframe so as to output
#     2023-Aug-25: move method to interpolate data to EmeraldMath.jl
#
#######################################################################################################################################################################################################
"""

    weather_driver(wd_tag::String, gmdict::Dict{String,Any}; appending::Bool = false)

Prepare weather driver dataframe to feed SPAC, given
- `wd_tag` Weather driver version tag
- `gmdict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add

"""
function weather_driver(wd_tag::String, gmdict::Dict{String,Any}; appending::Bool = false)
    nc_wd = weather_driver_file(wd_tag, gmdict; appending = appending)[1];
    df_wd = read_nc(nc_wd);

    # interpolate the data to a new resolution
    df_wd[!,"CO2"        ] .= interpolate_data(gmdict["CO2"], gmdict["YEAR"]; out_reso = "1H");
    df_wd[!,"CHLOROPHYLL"] .= interpolate_data(gmdict["CHLOROPHYLL"], gmdict["YEAR"]; out_reso = "1H");
    df_wd[!,"CI"         ] .= interpolate_data(gmdict["CLUMPING"], gmdict["YEAR"]; out_reso = "1H");
    df_wd[!,"LAI"        ] .= interpolate_data(gmdict["LAI"], gmdict["YEAR"]; out_reso = "1H");
    df_wd[!,"VCMAX25"    ] .= interpolate_data(gmdict["VCMAX25"], gmdict["YEAR"]; out_reso = "1H");

    # add the fields to store outputs
    for label in DF_VARIABLES
        df_wd[!,label] .= 0.0;
    end;
    for label in DF_SIMULATIONS
        df_wd[!,label] .= NaN;
    end;

    return df_wd
end;

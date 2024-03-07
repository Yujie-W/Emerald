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

    weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false)

Prepare weather driver dataframe to feed SPAC, given
- `wd_tag` Weather driver version tag
- `gm_dict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add

"""
function weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false)
    nc_wd = weather_driver_file(wd_tag, gm_dict; appending = appending)[1];
    df_wd = read_nc(nc_wd);

    # interpolate the data to a new resolution
    df_wd[!,"CO2"    ] .= interpolate_data(gm_dict["CO2"], gm_dict["YEAR"]; out_reso = "1H");
    df_wd[!,"CHL"    ] .= interpolate_data(gm_dict["CHLOROPHYLL"], gm_dict["YEAR"]; out_reso = "1H");
    df_wd[!,"CI"     ] .= interpolate_data(gm_dict["CLUMPING"], gm_dict["YEAR"]; out_reso = "1H");
    df_wd[!,"LAI"    ] .= interpolate_data(gm_dict["LAI"], gm_dict["YEAR"]; out_reso = "1H");
    df_wd[!,"VCMAX25"] .= interpolate_data(gm_dict["VCMAX25"], gm_dict["YEAR"]; out_reso = "1H");

    # add the fields to store outputs
    for label in DF_VARIABLES
        df_wd[!,label] .= 0.0;
    end;
    for label in DF_SIMULATIONS
        df_wd[!,label] .= NaN;
    end;

    return df_wd
end;

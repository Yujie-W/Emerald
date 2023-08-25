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

    weather_driver(wd_tag::String, gmdict::Dict{String,Any}; appending::Bool = false, displaying::Bool = true)

Prepare weather driver dataframe to feed SPAC, given
- `wd_tag` Weather driver version tag
- `gmdict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add
- `displaying` If true, display information about the NetCDF file

"""
function weather_driver(wd_tag::String, gmdict::Dict{String,Any}; appending::Bool = false, displaying::Bool = true)
    _nc_wd = weather_driver_file(wd_tag, gmdict; appending = appending, displaying = displaying)[1];
    _df_wd = read_nc(_nc_wd);

    # interpolate the data to a new resolution
    _df_wd[!,"CO2"        ] .= interpolate_data(gmdict["CO2"], gmdict["YEAR"]; out_reso = "1H");
    _df_wd[!,"CHLOROPHYLL"] .= interpolate_data(gmdict["CHLOROPHYLL"], gmdict["YEAR"]; out_reso = "1H");
    _df_wd[!,"CI"         ] .= interpolate_data(gmdict["CLUMPING"], gmdict["YEAR"]; out_reso = "1H");
    _df_wd[!,"LAI"        ] .= interpolate_data(gmdict["LAI"], gmdict["YEAR"]; out_reso = "1H");
    _df_wd[!,"VCMAX25"    ] .= interpolate_data(gmdict["VCMAX25"], gmdict["YEAR"]; out_reso = "1H");

    # add the fields to store outputs
    for _label in DF_VARIABLES
        _df_wd[!,_label] .= 0.0;
    end;
    for _label in DF_SIMULATIONS
        _df_wd[!,_label] .= NaN;
    end;

    return _df_wd
end

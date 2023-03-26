#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: move function from ClimaLand-0.2
#
#######################################################################################################################################################################################################
"""

    weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false, displaying::Bool = true)

Prepare weather driver dataframe to feed SPAC, given
- `wd_tag` Weather driver version tag
- `gm_dict` Dictionary that store grid information
- `appending` If true, always check whether there are new fields to add
- `displaying` If true, display information about the NetCDF file

"""
function weather_driver(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false, displaying::Bool = true)
    _nc_wd = weather_driver_file(wd_tag, gm_dict; appending = appending, displaying = displaying)[1];
    _df_wd = read_nc(_nc_wd);

    #
    # extropolate the time series based on input variable dimensions, dimensions must be within supported settings
    #     1. extropolate the data to 1D resolution
    #     2. extropolate the data to 1H resolution
    #
    _year = gm_dict["YEAR"];
    _days = isleapyear(_year) ? 366 : 365;
    @inline nt_to_1h(label::String) = (
        _dat_in = gm_dict[label];
        @assert length(_dat_in) in [366, 365, 53, 52, 46, 12, 1] "Dataset length not supported";

        if length(_dat_in) == 1
            _dat_1d = repeat([_dat_in;]; inner = _days);
        elseif length(_dat_in) == 12
            _dat_1d = [([repeat(_dat_in[_m:_m], month_days(_year, _m)) for _m in 1:12]...)...]
        elseif length(_dat_in) == 46
            _dat_1d = repeat(_dat_in; inner = 8)[1:_days]
        elseif length(_dat_in) in [52,53]
            _dat_1d = repeat([_dat_in;_dat_in[end]]; inner = 7)[1:_days]
        elseif length(_dat_in) in [365,366]
            _dat_1d = [_dat_in;_dat_in[end]][1:_days]
        end;

        return repeat(_dat_1d; inner = 24)
    );
    _df_wd[!,"CO2"        ] .= nt_to_1h("CO2");
    _df_wd[!,"CHLOROPHYLL"] .= nt_to_1h("CHLOROPHYLL");
    _df_wd[!,"CI"         ] .= nt_to_1h("CLUMPING");
    _df_wd[!,"LAI"        ] .= nt_to_1h("LAI");
    _df_wd[!,"VCMAX25"    ] .= nt_to_1h("VCMAX25");

    # add the fields to store outputs
    for _label in DF_VARIABLES
        _df_wd[!,_label] .= 0.0;
    end;

    return _df_wd
end

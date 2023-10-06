module Data

using Dates: isleapyear

using ..EmeraldUtility.Time: month_days


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Aug-25: add function (moved from EmeraldEarth.jl)
#
#######################################################################################################################################################################################################
"""

    interpolate_data!(data::Union{FT, Vector{FT}}) where {FT}

Gap fill the data linearly, given
- `data` Input data

"""
function interpolate_data!(data::Union{FT, Vector{FT}}) where {FT}
    if sum(.!isnan.(data)) in [0, length(data)]
        return nothing
    end;

    @inline find_last_number(vec_in::Vector{FT}, ind::Int) where {FT} = (
        _x = ind;
        _y = vec_in[ind];
        for _i in ind:-1:1
            if !isnan(vec_in[_i])
                _x = _i;
                _y = vec_in[_i];
                break;
            end;
        end;

        return _x, _y
    );

    @inline find_next_number(vec_in::Vector{FT}, ind::Int) where {FT} = (
        _x = ind;
        _y = vec_in[ind];
        for _i in ind:1:length(vec_in)
            if !isnan(vec_in[_i])
                _x = _i;
                _y = vec_in[_i];
                break;
            end;
        end;

        return _x, _y
    );

    @inline interpolate_data!(vec_in::Vector{FT}, ind::Int) where {FT} = (
        if isnan(vec_in[ind])
            (_x1,_y1) = find_last_number(vec_in, ind);
            (_x2,_y2) = find_next_number(vec_in, ind);
            vec_in[ind] = ((ind - _x1) * _y2 + (_x2 - ind) * _y1) / (_x2 - _x1);
        end;

        return nothing
    );

    _data_3x = [data; data; data];
    interpolate_data!.([_data_3x], (length(data)+1):(length(data)*2));
    data .= _data_3x[(length(data)+1):(length(data)*2)];

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Aug-25: add function (moved from EmeraldFrontier.jl)
#     2023-Aug-25: add support for single value number (not an array)
#
#######################################################################################################################################################################################################
"""

    interpolate_data(dat_in::Union{FT,Vector{FT}}, year::Int64; out_reso::String = "1H") where {FT}

Interpolate the data to 1H or 1D resolution, given
- `dat_in` Input data
- `year` Year of the input data
- `out_reso` Output temporal resolution

"""
function interpolate_data(dat_in::Union{FT,Vector{FT}}, year::Int64; out_reso::String = "1H") where {FT}
    @assert length(dat_in) in [366, 365, 53, 52, 46, 12, 1] "Dataset length not supported";

    _days = isleapyear(year) ? 366 : 365;

    if length(dat_in) == 1
        _dat_1d = repeat([dat_in;]; inner = _days);
    elseif length(dat_in) == 12
        _dat_1d = [([repeat(dat_in[_m:_m], month_days(year, _m)) for _m in 1:12]...)...];
    elseif length(dat_in) == 46
        _dat_1d = repeat(dat_in; inner = 8)[1:_days];
    elseif length(dat_in) in [52,53]
        _dat_1d = repeat([dat_in;dat_in[end]]; inner = 7)[1:_days];
    elseif length(dat_in) in [365,366]
        _dat_1d = [dat_in;dat_in[end]][1:_days];
    end;

    # return the data based on the output temporal resolution
    if out_reso == "1H"
        return repeat(_dat_1d; inner = 24)
    elseif out_reso == "1D"
        return _dat_1d
    end;
end;


end;

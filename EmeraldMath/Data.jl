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
function interpolate_data! end;

interpolate_data!(data::Union{FT, Vector{FT}}) where {FT} = (
    if sum(.!isnan.(data)) in [0, length(data)]
        return nothing
    end;

    data_3x = [data; data; data];
    interpolate_data!.([data_3x], (length(data)+1):(length(data)*2));
    data .= data_3x[(length(data)+1):(length(data)*2)];

    return nothing
);

interpolate_data!(vec_in::Vector{FT}, ind::Int) where {FT} = (
    if isnan(vec_in[ind])
        (xi,yi) = previous_number(vec_in, ind);
        (xj,yj) = next_number(vec_in, ind);
        vec_in[ind] = ((ind - xi) * yj + (xj - ind) * yi) / (xj - xi);
    end;

    return nothing
);

previous_number(vec_in::Vector{FT}, ind::Int) where {FT} = (
    xi = ind;
    yi = vec_in[ind];
    for i in ind:-1:1
        if !isnan(vec_in[i])
            xi = i;
            yi = vec_in[i];
            break;
        end;
    end;

    return xi, yi
);

next_number(vec_in::Vector{FT}, ind::Int) where {FT} = (
    xj = ind;
    yj = vec_in[ind];
    for j in ind:1:length(vec_in)
        if !isnan(vec_in[j])
            xj = j;
            yj = vec_in[j];
            break;
        end;
    end;

    return xj, yj
);


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

    nday = isleapyear(year) ? 366 : 365;

    if length(dat_in) == 1
        dat_1d = repeat([dat_in;]; inner = nday);
    elseif length(dat_in) == 12
        dat_1d = [([repeat(dat_in[_m:_m], month_days(year, _m)) for _m in 1:12]...)...];
    elseif length(dat_in) == 46
        dat_1d = repeat(dat_in; inner = 8)[1:nday];
    elseif length(dat_in) in [52,53]
        dat_1d = repeat([dat_in;dat_in[end]]; inner = 7)[1:nday];
    elseif length(dat_in) in [365,366]
        dat_1d = [dat_in;dat_in[end]][1:nday];
    end;

    # return the data based on the output temporal resolution
    if out_reso == "1H"
        return repeat(dat_1d; inner = 24)
    elseif out_reso == "1D"
        return dat_1d
    end;
end;


end;

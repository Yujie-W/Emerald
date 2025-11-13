module Data

using Dates: isleapyear

using ..Stats: nanmean
using ..EmeraldUtility.Time: month_days


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Aug-25: add function (moved from EmeraldEarth.jl)
#
#######################################################################################################################################################################################################
"""

    gapfill_data!(data::Union{FT, Vector{FT}}) where {FT}

Gap fill the data linearly, given
- `data` Input data

"""
function gapfill_data! end;

gapfill_data!(data::Union{FT, Vector{FT}}) where {FT} = (
    if sum(.!isnan.(data)) in [0, length(data)]
        return nothing
    end;

    data_3x = [data; data; data];
    gapfill_data!.([data_3x], (length(data)+1):(length(data)*2));
    data .= data_3x[(length(data)+1):(length(data)*2)];

    return nothing
);

gapfill_data!(vec_in::Vector{FT}, ind::Int) where {FT} = (
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
#     2024-Nov-13: move the method of read_spectrum as resample_data
#     2025-Nov-11: add supports for different output temporal resolutions (say 7D, 8D, and 1M)
#     2025-Nov-13: add supports for output temporal resolutions of 1Y
#
#######################################################################################################################################################################################################
"""

    resample_data(dat_in::Union{FT,Vector{FT}}, year::Int64; out_reso::String = "1H") where {FT}

Interpolate the data to 1H or 1D resolution, given
- `dat_in` Input data
- `year` Year of the input data
- `out_reso` Output temporal resolution

#

    resample_data(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT}

Linearly interpolate the data, given
- `x` Input x data
- `y` Input y data
- `target` Target x value

"""
function resample_data end;

resample_data(dat_in::Union{FT,Vector{FT}}, year::Int64; out_reso::String = "1H") where {FT} = (
    nday = isleapyear(year) ? 366 : 365;
    @assert length(dat_in) in [nday*24, nday, 53, 52, 46, 12, 1] "Dataset length not supported";
    @assert out_reso in ["1H", "1D", "7D", "8D", "1M", "1Y"] "Output temporal resolution not supported";

    dat_1d = if length(dat_in) == 1
        repeat([dat_in;]; inner = nday)
    elseif length(dat_in) == 12
        [([repeat(dat_in[_m:_m], month_days(year, _m)) for _m in 1:12]...)...]
    elseif length(dat_in) == 46
        repeat(dat_in; inner = 8)[1:nday]
    elseif length(dat_in) in [52,53]
        repeat([dat_in;dat_in[end]]; inner = 7)[1:nday]
    elseif length(dat_in) == nday
        dat_in
    elseif length(dat_in) == nday*24
        [nanmean(dat_in[((d8-1)*24+1):(d8*24)]) for d8 in 1:nday]
    end;

    # return the data based on the output temporal resolution
    return if out_reso == "1H"
        repeat(dat_1d; inner = 24)
    elseif out_reso == "1D"
        dat_1d
    elseif out_reso == "7D"
        [nanmean(dat_1d[((wk-1)*7+1):min(wk*7, nday)]) for wk in 1:53]
    elseif out_reso == "8D"
        [nanmean(dat_1d[((d8-1)*8+1):min(d8*8, nday)]) for d8 in 1:46]
    elseif out_reso == "1M"
        [nanmean(dat_1d[month_days(year, m; ranges = true)]) for m in 1:12]
    elseif out_reso == "1Y"
        nanmean(dat_1d)
    end;
);

resample_data(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT} = (
    @assert length(x) == length(y) "Dimensions of provided spectrum x and y must match!";
    @assert x[1] <= target <= x[end] "Target wavelength must be within the range provided spectum!";

    # iterate through the spectrum and find the index
    ind = 0;
    for i in 1:length(x)-1
        if x[i] <= target <= x[i+1]
            ind = i;
            break;
        end;
    end;

    return ((x[ind+1] - target) * y[ind] + (target - x[ind]) * y[ind+1]) / (x[ind+1] - x[ind])
);


end;

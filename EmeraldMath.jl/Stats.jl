module Stats

using Statistics: mean, median, std
using StatsBase: percentile


#######################################################################################################################################################################################################
#
# Changes to the functions
# General
#     2022-Oct-17: move functions outside of the folder
#     2023-Mar-10: move function to Emerald
#     2023-Mar-27: and method for nanmax, nanmin, etc to handle single number rather than an array
#     2023-Aug-04: add function nansum
#
#######################################################################################################################################################################################################
"""

    nanmax(x::Array)
    nanmax(x::Number)

Return the maximum of array ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmax end

nanmax(x::Array) = (
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : maximum( _x )
);

nanmax(x::Number) = isnan(x) ? NaN : x;


"""

    nanmean(x::Array)
    nanmean(x::Number)

Return the mean of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmean end

nanmean(x::Array) = (
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : mean( _x )
);

nanmean(x::Number) = isnan(x) ? NaN : x;


"""

    nanmedian(x::Array)
    nanmedian(x::Number)

Return the median of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmedian end

nanmedian(x::Array) = (
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : median( _x )
);

nanmedian(x::Number) = isnan(x) ? NaN : x;


"""

    nanmin(x::Array)
    nanmin(x::Number)

Return the maximum of array ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmin end
nanmin(x::Array) = (
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : minimum( _x )
);

nanmin(x::Number) = isnan(x) ? NaN : x;


"""

    nanpercentile(x::Array, p::Number)
    nanpercentile(x::Number, p::Number)

Return the percentile by excluding the NaN of given
- `x` Array of data
- `p` Percentile in `[%]`

"""
function nanpercentile end

nanpercentile(x::Array, p::Number) = (
    @assert 0 <= p <= 100;

    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : percentile( _x, p )
);

nanpercentile(x::Number, p::Number) = isnan(x) ? NaN : x;


"""

    nanstd(x::Array)
    nanstd(x::Number)

Return the std of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanstd end

nanstd(x::Array) = (
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : std( _x )
);

nanstd(x::Number) = isnan(x) ? NaN : x;


"""

    nansum(x::Array)
    nansum(x::Number)

Return the sum of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nansum end

nansum(x::Array) = (
    _x = filter(!isnan, x);

    return length(_x) == 0 ? 0 : sum( _x )
);

nansum(x::Number) = isnan(x) ? 0 : x;


#######################################################################################################################################################################################################
#
# Changes to the functions
# General
#     2022-Oct-17: move functions outside of the folder
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    mae(y::Array, pred::Array)

Return the mean absolute error by ommiting the NaN, given
- `y` Array of numbers, can be NaN
- `pred` Array of predictions, can be NaN

"""
function mae(y::Array, pred::Array)
    return nanmean( abs.(y .- pred) )
end


"""

    mape(y::Array, pred::Array)

Return the mean absolute percentage error by ommiting the NaN, given
- `y` Array of numbers, can be NaN
- `pred` Array of predictions, can be NaN

"""
function mape(y::Array, pred::Array)
    _mean = abs( nanmean(y) );
    _diff = abs.(y .- pred) ./ _mean .* 100;

    return nanmean( _diff )
end


"""

    mase(y::Array, pred::Array)

Return the mean absolute standardized error by ommiting the NaN, given
- `y` Array of numbers, can be NaN
- `pred` Array of predictions, can be NaN

"""
function mase(y::Array, pred::Array)
    _nstd = nanstd(y);
    _diff = abs.(y .- pred) ./ _nstd .* 100;

    return nanmean( _diff )
end


"""

    rmse(y::Array, pred::Array)

Return the root mean square error by ommiting the NaN, given
- `y` Array of numbers, can be NaN
- `pred` Array of predictions, can be NaN

"""
function rmse(y::Array, pred::Array)
    return sqrt( nanmean( (y .- pred) .^ 2 ) )
end


end

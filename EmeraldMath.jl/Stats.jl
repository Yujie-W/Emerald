module Stats

using Statistics: mean, median, std
using StatsBase: percentile


#######################################################################################################################################################################################################
#
# Changes to the functions
# General
#     2022-Oct-17: move functions outside of the folder
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    nanmax(x::Array)

Return the maximum of array ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmax(x::Array)
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : maximum( _x )
end


"""

    nanmean(x::Array)

Return the mean of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmean(x::Array)
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : mean( _x )
end


"""

    nanmedian(x::Array)

Return the median of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmedian(x::Array)
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : median( _x )
end


"""

    nanmin(x::Array)

Return the maximum of array ommiting the NaN, given
- `x` Array of numbers, can be NaN

"""
function nanmin(x::Array)
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : minimum( _x )
end


"""

    nanpercentile(x::Array, p::Number)

Return the percentile by excluding the NaN of given
- `x` Array of data
- `p` Percentile in `[%]`

"""
function nanpercentile(x::Array, p::Number)
    @assert 0 <= p <= 100

    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : percentile( _x, p )
end


"""

    nanstd(x::Array)

Return the std of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

```
"""
function nanstd(x::Array)
    _x = filter(!isnan, x);

    return length(_x) == 0 ? NaN : std( _x )
end


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

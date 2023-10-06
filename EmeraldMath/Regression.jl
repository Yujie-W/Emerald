module Regression

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using GLM: adjr², coef, coeftable, confint, lm, predict


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Nov-14: move function from PlotPlants
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save linear regression results

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LinearRegressionResult
    "Regression results"
    LM

    "Coeffients for all X"
    COEF::Vector = coef(LM)
    "Confidence intervals of the coefficients"
    CONFINT::Matrix = confint(LM)
    "P values of each coefficient"
    P::Vector = coeftable(LM).cols[4]
    "Adjust R²"
    R²::Number = adjr²(LM)
    "Predicted Y with confidence interval"
    XY::DataFrame = (
        _df = DataFrame();
        _pred = predict(LM, LM.pp.X; interval=:confidence);
        for i in axes(LM.pp.X,2)
            _df[!,"X$(i)"] = LM.pp.X[:,i];
        end;
        _df[!,"Y"     ] = LM.rr.y;
        _df[!,"predY" ] = _pred.prediction;
        _df[!,"lowerY"] = _pred.lower[:,1];
        _df[!,"upperY"] = _pred.upper[:,1];
        _df)
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-14: move function from PlotPlants
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    linear_regress(xs::Tuple, y::Vector)

Return the linear regression results, given
- `xs` A tuple of x variables
- `y` A vector of y

---
# Example
```julia
x1 = rand(5);
x2 = rand(5); x2[2] = NaN;
yy = rand(5);
lr1 = linear_regress((x1,), yy);
lr2 = linear_regress((x1,1), yy);
lr3 = linear_regress((x1,x2), yy);
lr4 = linear_regress((x1,x2,1), yy);
```

"""
function linear_regress(xs::Tuple, y::Vector)
    # make sure that the vectors match in dimensions
    for i in eachindex(xs)
        @assert (length(xs[i]) == 1) || (length(xs[i]) == length(y) > 1) "X and Y must have the same dimension";
    end;

    # create a mask to mask out NaN values
    _mask = .!isnan.(y);
    for i in eachindex(xs)
        if length(xs[i]) > 1
            _mask .= _mask .&& .!isnan.(xs[i]);
        end;
    end;

    # copy xs to _matx
    _matx = zeros(sum(_mask), length(xs));
    _vecy = zeros(sum(_mask));
    for i in eachindex(xs)
        if length(xs[i]) > 1
            _matx[:,i] .= xs[i][_mask];
        else
            _matx[:,i] .= xs[i];
        end;
    end;
    _vecy .= y[_mask];

    # run the fitting using GLM
    return LinearRegressionResult(LM = lm(_matx, _vecy))
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-14: move function from PlotPlants
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    test_slope(xs::Tuple, y::Vector; slope::Number = 0)

Return the P value of whether linear regression is same as the provided slope, given
- `xs` A tuple of x variables
- `y` A vector of y
- `slope` Target slope

---
# Example
```julia
x1 = rand(5);
yy = rand(5);
p1 = test_slope((x1,), yy; slope = 0);
p2 = test_slope((x1,1), yy; slope = 1);
```

"""
function test_slope(xs::Tuple, y::Vector; slope::Number = 0)
    # make sure that the vectors match in dimensions
    for i in eachindex(xs)
        @assert (length(xs[i]) == 1) || (length(xs[i]) == length(y) > 1) "X and Y must have the same dimension";
    end;
    @assert length(xs[1]) > 1 "First X cannot be 1";

    # create a mask to mask out NaN values
    _mask = .!isnan.(y);
    for i in eachindex(xs)
        if length(xs[i]) > 1
            _mask .= _mask .&& .!isnan.(xs[i]);
        end;
    end;

    # copy xs to _matx
    _matx = zeros(sum(_mask), length(xs));
    _vecy = zeros(sum(_mask));
    for i in eachindex(xs)
        if length(xs[i]) > 1
            _matx[:,i] .= xs[i][_mask];
        else
            _matx[:,i] .= xs[i];
        end;
    end;
    _vecy .= y[_mask] .- slope .* xs[1][_mask];

    # run the fitting using GLM
    _lr = lm(_matx, _vecy);

    return coeftable(_lr).cols[4][1]
end


end # module

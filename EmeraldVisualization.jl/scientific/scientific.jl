#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-14: refactor the function
#     2022-Nov-14: use EmeraldRegression function (and remove option intercept)
#
#######################################################################################################################################################################################################
"""

    plot_density!(ax::PyObject, xs::Vector, ys::Vector; cmap::String = "viridis", dmax::Union{Nothing, Number} = nothing, markersize::Number = 8)

Plot the scatter plot colored by the density of dots, given
- `ax` PyPlot Axis to plot on
- `xs` Vector of X
- `ys` Vector of Y
- `cmap` Color map scheme
- `dmax` Maximum density in the color scheme
- `markersize` Marker size (use the square value when plotting)

"""
function plot_density!(ax::PyObject, xs::Vector, ys::Vector; cmap::String = "viridis", dmax::Union{Nothing, Number} = nothing, markersize::Number = 8)
    # prepare dataframe
    _mask = (.!isnan.(xs)) .&& (.!isnan.(ys));
    _varx = xs[_mask];
    _vary = ys[_mask];
    _varc = similar(_varx);

    _ik = kde((_varx, _vary));
    if length(_varx) > 1000
        @showprogress for i in eachindex(_varx)
            _varc[i] = pdf(_ik, _varx[i], _vary[i]);
        end
    else
        _varc .= pdf.([_ik], _varx, _vary);
    end
    df = DataFrame(X = _varx, Y = _vary, C = _varc);
    sort!(df, [:C]);

    # set a max density limit for the data
    _c = isnothing(dmax) ? df.C : min.(df.C, dmax);
    ax.scatter(df.X, df.Y, s = markersize^2, c = _c, cmap = cmap);

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jan-19: refactor the function to plot_hexbin (because there is returned plot)
#
#######################################################################################################################################################################################################
"""

    plot_hexbin(ax::PyObject, xs::Array, ys::Array; cmap::String = "Greys", logbins::Bool = false, gridsize::Number = 25)
    plot_hexbin(ax::PyObject, xs::Array, ys::Array, xlims::Vector, ylims::Vector; cmap::String = "Greys", logbins::Bool = false, gridsize::Number = 25)

Return the density plot, given
- `ax` Axis to plot on
- `xs` Array of X
- `ys` Array of Y
- `cmap` Optional. Color map scheme
- `logbins` Optional. If true, use log(count) to color the bins
- `gridsize` Number of bins on both directions
- `xlim` Limits of x axis. Used to make plot region equal among subplots
- `ylim` Limits of y axis. Used to make plot region equal among subplots

"""
function plot_hexbin end

plot_hexbin(ax::PyObject, xs::Array, ys::Array; cmap::String = "Greys", logbins::Bool = false, gridsize::Number = 25) = (
    if logbins
        return ax.hexbin(xs, ys; cmap=cmap, bins="log", gridsize=gridsize)
    end;

    return ax.hexbin(xs, ys; cmap=cmap, gridsize=gridsize)
);

plot_hexbin(ax::PyObject, xs::Array, ys::Array, xlims::Vector, ylims::Vector; cmap::String = "Greys", logbins::Bool = false, gridsize::Number = 25) = (
    newxs = [collect(xs); xlims[1] - (xlims[2] - xlims[1]) / 10; xlims[2] + (xlims[2] - xlims[1]) / 10];
    newys = [collect(ys); ylims[1] - (ylims[2] - ylims[1]) / 10; ylims[2] + (ylims[2] - ylims[1]) / 10];

    return plot_hexbin(ax, newxs, newys; cmap=cmap, logbins=logbins, gridsize=gridsize)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-14: refactor the function
#     2022-Nov-14: use EmeraldRegression function (and remove option intercept)
#
#######################################################################################################################################################################################################
"""

    plot_line_regress!(ax::PyObject, xs::Tuple, ys::Vector; linestyle::String = "-", interval::Bool = false, color::String = "red", alpha::Number = 0.3)

Plor linear regression and confidence interval on the axis, given
- `ax` Given axis
- `xs` Tuple of x
- `ys` Vector of y
- `linestyle` Optional. Line style for the regression curve ("-" by default)
- `interval` Optional: if true, plot the confidence interval of fitted y
- `color` Color the fitted curve
- `alpha` Transparency of the confidence interval (same color as curve)

"""
function plot_line_regress!(ax::PyObject, xs::Tuple, ys::Vector; linestyle::String = "-", interval::Bool = false, color::String = "red", alpha::Number = 0.3)
    # make linear regression
    lr = linear_regress(xs, ys);
    df = sort(lr.XY);

    # plot the fittings
    ax.plot(df.X1, df.predY; linestyle = linestyle, color = color);
    if interval
        ax.fill_between(df.X1, df.lowerY, df.upperY; facecolor = color, alpha = alpha);
    end

    return nothing
end

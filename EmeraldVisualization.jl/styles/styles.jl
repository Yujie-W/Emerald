#######################################################################################################################################################################################################
#
# Changes to these types and structs
# General
#     2023-Jan-20: add new function
#
#######################################################################################################################################################################################################
"""

Abstract type of subplot styles used for visualization for Plots.jl

$(TYPEDEF)

"""
abstract type AbstractSubplotStyle end


"""

Style for plotting global maps

$(TYPEDEF)

"""
struct GlobalMapStyle <: AbstractSubplotStyle end


"""

Style for plotting regular figures with legends

$(TYPEDEF)

"""
struct RegularPlotStyle <: AbstractSubplotStyle end


"""

Style for plotting regular figures without legends

$(TYPEDEF)

"""
struct RegularPlotStyleNoLegend <: AbstractSubplotStyle end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jan-20: add new function
#
#######################################################################################################################################################################################################
"""

    set_style!(fig::Plot, style::GlobalMapStyle; latitude_360::Bool = false)
    set_style!(panel::Subplot, style::GlobalMapStyle; latitude_360::Bool = false)
    set_style!(fig::Plot, style::Union{RegularPlotStyle, RegularPlotStyleNoLegend})
    set_style!(panel::Subplot, style::Union{RegularPlotStyle, RegularPlotStyleNoLegend})

Set the figure or subplot style, given
- `fig` Plots.Plot type Figure
- `style` Subplot style struct
- `panel` Subplot panel
- `latitude_360` Option for GlobalMapStyle. If true, use 0-360 for latitude

"""
function set_style! end

set_style!(fig::Plot, style::GlobalMapStyle; latitude_360::Bool = false) = set_style!.(fig.subplots, [style]; latitude_360 = latitude_360);

set_style!(panel::Subplot, style::GlobalMapStyle; latitude_360::Bool = false) = (
    panel[:framestyle] = :box;
    panel[:aspect_ratio] = :equal;

    # set up the lims and ticks
    (_xmin,_xmax) = latitude_360 ? (0,360) : (-180,180);
    _xticklabels = latitude_360 ? String[string(_lon) for _lon in _xmin:60:_xmax] : String["180W", "120W", "60W", "0", "60E", "120E", "180E"];
    panel[:xaxis][:lims] = (_xmin,_xmax);
    panel[:yaxis][:lims] = (-90,90);
    panel[:xaxis][:ticks] = (collect(_xmin:60:_xmax), _xticklabels);
    panel[:yaxis][:ticks] = (collect(-90:30:90), String["90S", "60S", "30S", "0", "30N", "60N", "90N"]);

    return nothing
);

set_style!(fig::Plot, style::Union{RegularPlotStyle, RegularPlotStyleNoLegend}) = set_style!.(fig.subplots, [style]);

set_style!(panel::Subplot, style::Union{RegularPlotStyle, RegularPlotStyleNoLegend}) = (
    panel[:framestyle] = :box;
    panel[:legend_position] = (style isa RegularPlotStyle ? :best : :none);
    panel[:legend_font_pointsize] = 15;

    for _axis in [:xaxis, :yaxis, :zaxis]
        panel[_axis][:gridlinewidth] = 1;
        panel[_axis][:gridstyle] = :dash;
        panel[_axis][:tickfontsize] = 15;
    end;

    return nothing
);

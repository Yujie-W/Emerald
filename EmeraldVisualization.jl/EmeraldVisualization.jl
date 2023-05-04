module EmeraldVisualization

import Plots

using DataFrames: DataFrame, sort!
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using KernelDensity: kde, pdf
using Plots: Animation, Plot, Subplot, gif, heatmap, @animate
using ProgressMeter: @showprogress
using PyCall: PyObject, pyimport
using PyPlot: Figure, figure, rc

using ..EmeraldIO.Netcdf: read_nc
using ..EmeraldMath.Regression: linear_regress


PATCHES = pyimport("matplotlib.patches");


include("animation/animation.jl");
include("animation/netcdf.jl");

include("scientific/canvas.jl");
include("scientific/decoration.jl");
include("scientific/latex.jl");
include("scientific/scientific.jl");
include("scientific/shapes.jl");

include("styles/styles.jl");


end # module

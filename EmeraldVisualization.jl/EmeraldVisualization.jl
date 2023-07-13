module EmeraldVisualization

import Plots

using DataFrames: DataFrame, sort!
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using KernelDensity: kde, pdf
using Plots: Animation, Plot, Subplot, gif, heatmap, savefig, @animate
using ProgressMeter: @showprogress
using PyCall: PyObject, pyimport
using PyPlot: Figure, figure, rc

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: read_nc

using ..EmeraldMath.Regression: linear_regress
using ..EmeraldMath.Stats: nanmin, nanmax


PATCHES = pyimport("matplotlib.patches");

LAND_MASK = nothing;


include("animation/animation.jl");
include("animation/data.jl");
include("animation/netcdf.jl");

include("data/land_mask.jl");

include("scientific/canvas.jl");
include("scientific/decoration.jl");
include("scientific/latex.jl");
include("scientific/scientific.jl");
include("scientific/shapes.jl");

include("styles/styles.jl");


end # module

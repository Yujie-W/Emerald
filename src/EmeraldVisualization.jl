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


include("../EmeraldVisualization.jl/animation/animation.jl");
include("../EmeraldVisualization.jl/animation/data.jl");
include("../EmeraldVisualization.jl/animation/netcdf.jl");

include("../EmeraldVisualization.jl/data/land_mask.jl");

include("../EmeraldVisualization.jl/scientific/canvas.jl");
include("../EmeraldVisualization.jl/scientific/decoration.jl");
include("../EmeraldVisualization.jl/scientific/latex.jl");
include("../EmeraldVisualization.jl/scientific/scientific.jl");
include("../EmeraldVisualization.jl/scientific/shapes.jl");

include("../EmeraldVisualization.jl/styles/styles.jl");


end # module

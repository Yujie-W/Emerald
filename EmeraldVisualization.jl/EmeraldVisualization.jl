module EmeraldVisualization

import Plots

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using Plots: Animation, Plot, Subplot, gif, heatmap, @animate

using ..EmeraldIO.Netcdf: read_nc


include("animation/animation.jl");
include("animation/netcdf.jl");

include("styles/styles.jl");


end # module

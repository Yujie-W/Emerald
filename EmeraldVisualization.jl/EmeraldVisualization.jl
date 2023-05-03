module EmeraldVisualization

using Plots

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using Plots: Plot, Subplot

using ..EmeraldIO.Netcdf: read_nc


include("animation.jl");
include("netcdf.jl");
include("styles.jl");


end # module

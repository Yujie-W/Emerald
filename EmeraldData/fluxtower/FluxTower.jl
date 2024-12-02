module FluxTower

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using NetcdfIO: read_nc, save_nc!

using ..EmeraldIO.Folders: AMERIFLUX_DATA, AMERIFLUX_REPROCESSED, FLUXNET2015_DATA, FLUXNET2015_REPROCESSED
using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmean
using ..EmeraldPhysics.Constant: K_STEFAN, T₀


# include the general data struct and constructors
include("towers.jl");
include("fluxtowerdata.jl");
include("process.jl");
include("query.jl");


end # module

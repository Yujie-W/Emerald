module FluxTower

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using NetcdfIO: save_nc!

using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmean
using ..EmeraldPhysics.Constant: K_STEFAN, T₀


# file locations for AmeriFlux and FluxNet data
AMERIFLUX_FOLDER = "/home/wyujie/DATASERVER/field/AmeriFlux";
AMERIFLUX_DATA = "$(AMERIFLUX_FOLDER)/original";
AMERIFLUX_REPROCESSED = "$(AMERIFLUX_FOLDER)/reprocessed";
FLUXNET_FOLDER = "/home/wyujie/DATASERVER/field/Fluxnet2015";
FLUXNET_DATA = "$(FLUXNET_FOLDER)/data";
FLUXNET_REPROCESSED = "$(FLUXNET_FOLDER)/reprocessed";


# include the general data struct and constructors
include("fluxtower/towers.jl");
include("fluxtower/fluxtowerdata.jl");
include("fluxtower/process.jl");


end # module

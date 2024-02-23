module EmeraldData


using ..EmeraldIO
using ..EmeraldLand
using ..EmeraldMath
using ..EmeraldPhysics
using ..EmeraldUtility


include("ERA5.jl");
include("FluxTower.jl");
include("GlobalDatasets.jl");
include("GriddingMachineData.jl");


end; # module

module EmeraldLand

using ..EmeraldIO
using ..EmeraldMath
using ..EmeraldPhysics


include("Namespace/Namespace.jl");

include("EarthGeometry.jl");
include("PhysicalChemistry.jl");

include("LeafOptics/LeafOptics.jl");
include("CanopyOptics/CanopyOptics.jl");
include("Photosynthesis/Photosynthesis.jl");
include("SoilHydraulics/SoilHydraulics.jl");
include("PlantHydraulics/PlantHydraulics.jl");
include("StomatalModels/StomatalModels.jl");
include("SPAC/SPAC.jl");


end

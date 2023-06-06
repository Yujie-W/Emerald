module EmeraldLand

using ..EmeraldIO
using ..EmeraldMath


include("Constant.jl");
include("Namespace.jl");

include("EarthGeometry.jl");
include("Optics.jl");
include("PhysicalChemistry.jl");

include("LeafOptics.jl");
include("CanopyOptics.jl");
include("Photosynthesis.jl");
include("SoilHydraulics.jl");
include("PlantHydraulics.jl");
include("StomatalModels.jl");
include("SPAC.jl");


end

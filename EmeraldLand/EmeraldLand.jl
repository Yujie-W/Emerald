module EmeraldLand

using ..EmeraldIO
using ..EmeraldMath


include("Constant.jl");
include("namespace/Namespace.jl");

include("EarthGeometry.jl");
include("Optics.jl");
include("PhysicalChemistry.jl");

include("optics/LeafOptics.jl");
include("radiation/CanopyOptics.jl");
include("photosynthesis/Photosynthesis.jl");
include("soil/SoilHydraulics.jl");
include("PlantHydraulics/PlantHydraulics.jl");
include("stomata/StomatalModels.jl");
include("spac/SPAC.jl");


end

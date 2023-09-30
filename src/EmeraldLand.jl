module EmeraldLand

using ..EmeraldIO
using ..EmeraldMath


include("../EmeraldLand.jl/Constant.jl");
include("../EmeraldLand.jl/namespace/Namespace.jl");

include("../EmeraldLand.jl/EarthGeometry.jl");
include("../EmeraldLand.jl/Optics.jl");
include("../EmeraldLand.jl/PhysicalChemistry.jl");

include("../EmeraldLand.jl/LeafOptics.jl");
include("../EmeraldLand.jl/CanopyOptics.jl");
include("../EmeraldLand.jl/Photosynthesis.jl");
include("../EmeraldLand.jl/SoilHydraulics.jl");
include("../EmeraldLand.jl/hydraulics/PlantHydraulics.jl");
include("../EmeraldLand.jl/StomatalModels.jl");
include("../EmeraldLand.jl/spac/SPAC.jl");


end

module EmeraldLand

using ..EmeraldIO
using ..EmeraldMath


include("../EmeraldLand.jl/Constant.jl");
include("../EmeraldLand.jl/Namespace.jl");

include("../EmeraldLand.jl/EarthGeometry.jl");
include("../EmeraldLand.jl/Optics.jl");
include("../EmeraldLand.jl/PhysicalChemistry.jl");

include("../EmeraldLand.jl/LeafOptics.jl");
include("../EmeraldLand.jl/CanopyOptics.jl");
include("../EmeraldLand.jl/Photosynthesis.jl");
include("../EmeraldLand.jl/SoilHydraulics.jl");
include("../EmeraldLand.jl/PlantHydraulics.jl");
include("../EmeraldLand.jl/StomatalModels.jl");
include("../EmeraldLand.jl/SPAC.jl");


end
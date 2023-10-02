module EmeraldLand

using ..EmeraldIO
using ..EmeraldMath
using ..EmeraldPhysics


# no internal inter-dependencies
include("Namespace/Namespace.jl");

# depends on Namespace
include("CanopyOptics/CanopyOptics.jl");
include("LeafOptics/LeafOptics.jl");
include("Photosynthesis/Photosynthesis.jl");
include("PhysicalChemistry.jl");

# depends on PhysicalChemistry
include("SoilHydraulics/SoilHydraulics.jl");

# depends on SoilHydraulics
include("PlantHydraulics/PlantHydraulics.jl");

# depends on PlantHydraulics
include("EnergyBudget/EnergyBudget.jl");
include("StomatalModels/StomatalModels.jl");

# depends on StomatalModels
include("SPAC/SPAC.jl");


end

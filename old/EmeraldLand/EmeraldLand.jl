module EmeraldLand

using ..EmeraldIO


# no internal inter-dependencies
include("Namespace/Namespace.jl");

# depends on Namespace
include("LeafOptics/LeafOptics.jl");
include("Photosynthesis/Photosynthesis.jl");
include("SoilHydraulics/SoilHydraulics.jl");

# depends on LeafOptics
include("CanopyOptics/CanopyOptics.jl");

# depends on SoilHydraulics
include("PlantHydraulics/PlantHydraulics.jl");

# depends on PlantHydraulics
include("EnergyBudget/EnergyBudget.jl");
include("StomatalModels/StomatalModels.jl");

# depends on StomatalModels
include("SPAC/SPAC.jl");


end;

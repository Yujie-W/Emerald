module Emerald

using Revise


# requires DataCenter to read the input data
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

# depends on SPAC
include("Land/Land.jl");


end # module Emerald

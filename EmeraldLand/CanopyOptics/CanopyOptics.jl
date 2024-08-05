module CanopyOptics

using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using SpecialFunctions: beta_inc
using Statistics: mean

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using ..EmeraldPhysics.Constant: K_STEFAN
using ..EmeraldPhysics.Optics: energy!, photon, photon!

using ..Namespace: BetaLIDF, VerhoefLIDF
using ..Namespace: CanopyStructureTrait, CanopyStructureTDAuxil, SensorGeometryState, SensorGeometrySDAuxil, SunGeometryState, SunGeometrySDAuxil
using ..Namespace: SoilAlbedoBroadbandCLIMA, SoilAlbedoBroadbandCLM, SoilAlbedoHyperspectralCLIMA, SoilAlbedoHyperspectralCLM, SoilAlbedoPrescribe
using ..Namespace: SoilLayer, SoilBulk
using ..Namespace: MultiLayerCanopy
using ..Namespace: CanopyLayer, Leaf
using ..Namespace: BulkSPAC, SPACConfiguration


# functions related to canopy geometry
include("geometry/direction.jl");
include("geometry/extinction.jl");
include("geometry/inclination.jl");
include("geometry/sensor.jl");
include("geometry/soil_albedo.jl");
include("geometry/structure.jl");
include("geometry/sun.jl");

# function related to canopy radiation
include("radiation/longwave.jl");
include("radiation/shortwave.jl");

include("radiation/fluorescence.jl");
include("radiation/reflection.jl");

include("radiation/pipeline.jl");


end; # module

module CanopyOptics

using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using SpecialFunctions: beta_inc
using Statistics: mean

using PkgUtility.MathTools: ReduceStepMethodND, SolutionToleranceND, find_peak
using PkgUtility.UniversalConstants: K_STEFAN
using PkgUtility.UniversalConstants: energy_to_photon, energy_to_photon!, photon_to_energy!

using ..Namespace: BetaLIDF, VerhoefLIDF
using ..Namespace: ClumpingIndex
using ..Namespace: CanopyStructureTrait, CanopyStructureTDAuxil, SensorGeometryState, SensorGeometrySDAuxil, SunGeometryState, SunGeometrySDAuxil
using ..Namespace: SoilAlbedoBroadbandCLIMA, SoilAlbedoBroadbandCLM, SoilAlbedoHyperspectralAsh, SoilAlbedoHyperspectralCLIMA, SoilAlbedoHyperspectralCLM, SoilAlbedoPrescribe
using ..Namespace: SoilLayer, SoilBulk
using ..Namespace: MultiLayerCanopy
using ..Namespace: CanopyLayer, Leaf
using ..Namespace: BulkSPAC, SPACConfig

using ..LeafOptics: layer_2_ρ, layer_2_τ, leaf_ρ, leaf_τ


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

module CanopyOptics

using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using SpecialFunctions: beta_inc
using Statistics: mean

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using ..EmeraldPhysics.Constant: K_STEFAN
using ..EmeraldPhysics.Optics: energy!, photon, photon!

using ..Namespace: MultiLayerCanopy, ReferenceSpectra, ShortwaveRadiation
using ..Namespace: BetaLIDF, Leaf, SoilBulk, SoilLayer, VerhoefLIDF
using ..Namespace: MultiLayerSPAC, SPACConfiguration


# functions related to canopy geometry
include("geometry/extinction.jl");
include("geometry/inclination.jl");
include("geometry/sensor.jl");
include("geometry/sun.jl");
include("geometry/structure.jl");

include("fluorescence.jl");
include("geometry.jl");
include("radiation.jl");
include("remote_sensing.jl");
include("soil.jl");


end; # module

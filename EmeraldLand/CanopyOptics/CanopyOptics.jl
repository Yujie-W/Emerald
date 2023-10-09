module CanopyOptics

using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using SpecialFunctions: beta_inc
using Statistics: mean

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using ..EmeraldPhysics.Constant: K_STEFAN
using ..EmeraldPhysics.Optics: energy!, photon, photon!

using ..Namespace: HyperspectralMLCanopy, ReferenceSpectra, ShortwaveRadiation
using ..Namespace: BetaLIDF, Leaf, SoilBulk, SoilLayer, SunSensorGeometry, VerhoefLIDF
using ..Namespace: MultiLayerSPAC, SPACConfiguration


include("clumping.jl");
include("coefficients.jl");
include("fluorescence.jl");
include("geometry.jl");
include("inclination.jl");
include("radiation.jl");
include("remote_sensing.jl");
include("soil.jl");


end # module

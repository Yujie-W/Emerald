module CanopyOptics

using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using SpecialFunctions: beta_inc
using Statistics: mean

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..Constant: K_STEFAN
using ..Namespace: BroadbandRadiation, BroadbandSLCanopy, BroadbandSoilAlbedo
using ..Namespace: HyperspectralMLCanopy, HyperspectralRadiation, HyperspectralSoilAlbedo, ReferenceSpectra
using ..Namespace: BetaLIDF, Leaves2D, Soil, SunSensorGeometry, VerhoefLIDF
using ..Namespace: MultiLayerSPAC, SPACConfiguration
using ..Optics: energy!, photon, photon!


include("clumping.jl");
include("coefficients.jl");
include("fluorescence.jl");
include("geometry.jl");
include("inclination.jl");
include("radiation.jl");
include("remote_sensing.jl");
include("soil.jl");


end # module

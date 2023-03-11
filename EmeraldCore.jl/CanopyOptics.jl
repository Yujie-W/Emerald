module CanopyOptics

using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using Statistics: mean

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..Constant: K_STEFAN
using ..Namespace: BroadbandRadiation, BroadbandSLCanopy, BroadbandSoilAlbedo
using ..Namespace: HyperspectralMLCanopy, HyperspectralRadiation, HyperspectralSoilAlbedo
using ..Namespace: Leaves1D, Leaves2D, Soil, SunSensorGeometry, VerhoefLIDF
using ..Namespace: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ..Optics: energy!, photon, photon!


include("radiation/clumping.jl"      )
include("radiation/coefficients.jl"  )
include("radiation/constants.jl"     )
include("radiation/fluorescence.jl"  )
include("radiation/geometry.jl"      )
include("radiation/inclination.jl"   )
include("radiation/radiation.jl"     )
include("radiation/remote_sensing.jl")
include("radiation/soil.jl"          )


end # module
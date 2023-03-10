module Solver

using DocStringExtensions: TYPEDEF, TYPEDFIELDS


include("solver/method.jl"   )
include("solver/tolerance.jl")

include("solver/find_peak.jl")
include("solver/find_zero.jl")


end

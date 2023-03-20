module Emerald

using Revise

export EmeraldCore
export EmeraldEarth
export EmeraldMath
export EmeraldIO
export EmeraldTest
export EmeraldUtility


# include the submodules
include("../EmeraldIO.jl/EmeraldIO.jl")
include("../EmeraldMath.jl/EmeraldMath.jl")
include("../EmeraldTest.jl/EmeraldTest.jl")
include("../EmeraldUtility.jl/EmeraldUtility.jl")
include("../EmeraldVisualization.jl/EmeraldVisualization.jl")

include("../EmeraldCore.jl/EmeraldCore.jl")
include("../EmeraldFrontier.jl/EmeraldFrontier.jl")
include("../EmeraldEarth.jl/EmeraldEarth.jl")


end # module

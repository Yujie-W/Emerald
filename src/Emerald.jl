module Emerald

export EmeraldCore
export EmeraldMath
export EmeraldIO
export EmeraldTest


# include the submodules
include("../EmeraldIO.jl/EmeraldIO.jl")
include("../EmeraldMath.jl/EmeraldMath.jl")
include("../EmeraldTest.jl/EmeraldTest.jl")
include("../EmeraldUtility.jl/EmeraldUtility.jl")
include("../EmeraldVisualization.jl/EmeraldVisualization.jl")

include("../EmeraldCore.jl/EmeraldCore.jl")


end # module

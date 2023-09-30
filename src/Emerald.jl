module Emerald

using Revise

export EmeraldData
export EmeraldEarth
export EmeraldFrontier
export EmeraldIO
export EmeraldLand
export EmeraldMath
export EmeraldTest
export EmeraldUtility


# include the submodules
include("../EmeraldIO.jl/EmeraldIO.jl");
include("../EmeraldTest.jl/EmeraldTest.jl");
include("../EmeraldUtility.jl/EmeraldUtility.jl");

include("../EmeraldMath.jl/EmeraldMath.jl");

include("../EmeraldLand.jl/EmeraldLand.jl");

include("../EmeraldData.jl/EmeraldData.jl");

include("../EmeraldEarth.jl/EmeraldEarth.jl");
include("../EmeraldFrontier.jl/EmeraldFrontier.jl");


end # module

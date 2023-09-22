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
include("EmeraldIO.jl");
include("EmeraldTest.jl");
include("EmeraldUtility.jl");

include("EmeraldMath.jl");

include("EmeraldLand.jl");

include("EmeraldData.jl");

include("EmeraldEarth.jl");
include("EmeraldFrontier.jl");


end # module

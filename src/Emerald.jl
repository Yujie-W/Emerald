module Emerald

using Revise

export EmeraldData
export EmeraldEarth
export EmeraldFrontier
export EmeraldIO
export EmeraldLand
export EmeraldMath
export EmeraldPhysics
export EmeraldTest
export EmeraldUtility


# include the submodules
include("../EmeraldIO/EmeraldIO.jl");
include("../EmeraldPhysics/EmeraldPhysics.jl");
include("../EmeraldTest/EmeraldTest.jl");
include("../EmeraldUtility/EmeraldUtility.jl");

include("../EmeraldMath/EmeraldMath.jl");

include("../EmeraldLand/EmeraldLand.jl");

include("../EmeraldData/EmeraldData.jl");

include("../EmeraldEarth/EmeraldEarth.jl");
include("../EmeraldFrontier/EmeraldFrontier.jl");


end # module

module Emerald

using Revise

export EmeraldData
export EmeraldEarth
export EmeraldFrontier
export EmeraldIO
export EmeraldLand


# include the submodules
include("../EmeraldIO/EmeraldIO.jl");

include("../EmeraldMath/EmeraldMath.jl");

include("../EmeraldLand/EmeraldLand.jl");
include("../EmeraldOcean/EmeraldOcean.jl");

include("../EmeraldData/EmeraldData.jl");

include("../EmeraldEarth/EmeraldEarth.jl");
include("../EmeraldFrontier/EmeraldFrontier.jl");


end; # module

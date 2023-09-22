using Emerald
using Test


@testset verbose = true "Emerald" begin
    #include("coverage/plant_hydraulics.jl");

    include("unittest.jl");

    #include("tutorial.jl");
end

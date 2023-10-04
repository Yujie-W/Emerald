using Emerald
using Test


@testset verbose = true "Emerald" begin
    include("unittest.jl");

    include("tutorial.jl");
end

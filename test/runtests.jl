using Emerald
using Test


@testset verbose = true "Emerald" begin
    include("core.jl");
    include("planthydraulics.jl");
end

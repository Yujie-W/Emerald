using Test


@testset verbose = true "Unit Test" begin
    include("unittests/canopyoptics.jl");
    include("unittests/leafoptics.jl");
    include("unittests/photosynthesis.jl");
    include("unittests/planthydraulics.jl");
    include("unittests/soilhydraulics.jl");
end;

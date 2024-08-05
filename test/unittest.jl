using Test


@testset verbose = true "Unit Test" begin
    include("unittests/namespace.jl");
    include("unittests/leafoptics.jl");
    include("unittests/canopyoptics.jl");
    include("unittests/photosynthesis.jl");
    include("unittests/soilhydraulics.jl");
    include("unittests/planthydraulics.jl");
    include("unittests/stomatalmodels.jl");
end;

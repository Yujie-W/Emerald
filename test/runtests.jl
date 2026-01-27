using Test


@testset verbose = true "Emerald" verbose = true begin
    @testset "Continous Integration" verbose = true begin
        include("Land/site-run.jl");
    end;

    @testset "ResearchTools" verbose = true begin
        include("ResearchTools/fit-aci.jl");
    end;

    @testset "Tutorials" verbose = true begin
        include("Tutorials/leaf-level-fluorescence.jl");
        include("Tutorials/leaf-level-photosynthesis.jl");
    end;
end;

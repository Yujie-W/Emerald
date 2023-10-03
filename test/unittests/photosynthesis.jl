using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS


@testset verbose = true "Plant Hydraulics Model" begin
    @testset "C3 Photosynthesis Model" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf{Float64}(config);
        leaf.NS.photosystem = NS.C3VJP{Float64}();
    end;
end;

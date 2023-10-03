using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS


@testset verbose = true "Plant Hydraulics Model" begin
    @testset "Temperature Dependencies" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            @test true;
        end;
    end;

    @testset "Electron Transport" begin
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            @test true;
        end;
    end;

end;

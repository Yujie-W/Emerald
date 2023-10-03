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
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            @test true;
        end;
    end;

    @testset "RubisCO Limited Rate (P Mode)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, 20.0);
            @test true;
        end;
    end;

    @testset "RubisCO Limited Rate (G Mode)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, air, 0.2);
            @test true;
        end;
    end;

    @testset "Light Limited Rate (P Mode)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.light_limited_rate!(ps);
            @test true;
        end;
    end;

    @testset "Light Limited Rate (G Mode)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.light_limited_rate!(ps, air, 0.2);
            @test true;
        end;
    end;

    @testset "Product Limited Rate (P Mode)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.product_limited_rate!(ps, 20.0);
            @test true;
        end;
    end;

    @testset "Product Limited Rate (G Mode)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.product_limited_rate!(ps, air, 0.2);
            @test true;
        end;
    end;

    @testset "Photosynthesis colimitation" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, 20.0);
            PS.light_limited_rate!(ps);
            PS.product_limited_rate!(ps, 20.0);
            PS.colimit_photosynthesis!(ps);
            @test true;
        end;
    end;

    @testset "Fluorescence coefficients" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, 20.0);
            PS.light_limited_rate!(ps);
            PS.product_limited_rate!(ps, 20.0);
            PS.colimit_photosynthesis!(ps);
            PS.photosystem_coefficients!(ps, 100.0);
            @test true;
        end;
    end;

end;

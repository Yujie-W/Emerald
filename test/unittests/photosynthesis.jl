using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS


@testset verbose = true "Photosynthesis Model" begin
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

    @testset "Photosynthesis only (for stomatal models)" begin
        air = NS.AirLayer{Float64}();
        for ps in [NS.C3VJP{Float64}(), NS.C4VJP{Float64}(), NS.C3Cyto{Float64}()]
            PS.photosynthesis_only!(ps, air, 0.2, 100.0, 298.15);
            @test true;
        end;
    end;

    @testset "Derivatives (for stomatal models)" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        @test PS.∂R∂T(leaf) > 0;

        for td in [NS.Arrhenius{Float64}(T_REF = 298.15, VAL_REF = NaN , ΔHA = 46390.0),
                   NS.ArrheniusPeak{Float64}(T_REF = 298.15, VAL_REF = NaN , ΔHA = 46390.0, ΔHD = 150650.0, ΔSV = 490.0),
                   NS.Q10{Float64}(Q_10 = 1.4, T_REF = 298.15, VAL_REF = 0.0140/8760),
                   NS.Q10Peak{Float64}(Q_10 = 1.4, T_REF = 298.15, VAL_REF = 0.0140/8760, ΔHD = 150650.0, ΔSV = 490.0)]
            @test PS.∂R∂T(td, 1.0, 298.15) > 0;
        end;
    end;

    @testset "Leaf Photosynthesis (Respiration only)" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = true);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(), 1.0; rd_only = true);
        @test true;
    end;

    @testset "Leaf Photosynthesis (core)" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar_sunlit .= 100.0;
        leaf.flux.auxil.ppar_shaded = 100.0;
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(), 1.0; rd_only = false);
        @test true;
    end;

    @testset "Leaf Photosynthesis (stomtal model)" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar_sunlit .= 100.0;
        leaf.flux.auxil.ppar_shaded = 100.0;

        leaf.flux.state.stomatal_model = NS.WangSM{Float64}();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(); rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(); rd_only = false);
        @test true;

        leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
        leaf.flux.state.stomatal_model.β.PARAM_Y = NS.BetaParameterG1();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(); rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(); rd_only = false);
        @test true;

        leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
        leaf.flux.state.stomatal_model.β.PARAM_Y = NS.BetaParameterVcmax();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(); rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(); rd_only = false);
        @test true;
    end;

    @testset "Plant Photosynthesis" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);

        PS.plant_photosynthesis!(spac, NS.GCO₂Mode());
        @test true;
        PS.plant_photosynthesis!(spac, NS.PCO₂Mode());
        @test true;

        spac = NS.MultiLayerSPAC(config);
        spac.CANOPY.lai = 0.0;
        PS.plant_photosynthesis!(spac, NS.GCO₂Mode());
        @test true;
        PS.plant_photosynthesis!(spac, NS.PCO₂Mode());
        @test true;

        spac = NS.MultiLayerSPAC(config);
        spac.ANGLES.sza = 90;
        PS.plant_photosynthesis!(spac, NS.GCO₂Mode());
        @test true;
        PS.plant_photosynthesis!(spac, NS.PCO₂Mode());
        @test true;
    end;
end;

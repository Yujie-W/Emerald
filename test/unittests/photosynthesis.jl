using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS
import Emerald.EmeraldLand.SPAC


@testset verbose = true "Photosynthesis.jl" begin
    @testset "Temperature dependencies" begin
        air = NS.AirLayer{Float64}();
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            @test true;
        end;
    end;

    @testset "Electron transport" begin
        air = NS.AirLayer{Float64}();
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            @test true;
        end;
    end;

    @testset "RubisCO limited rate" begin
        air = NS.AirLayer{Float64}();
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, 20.0);
            @test true;

            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, air, 0.2);
            @test true;
        end;
    end;

    @testset "Light limited rate" begin
        air = NS.AirLayer{Float64}();
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.light_limited_rate!(ps);
            @test true;

            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.light_limited_rate!(ps, air, 0.2);
            @test true;
        end;
    end;

    @testset "Product limited rate" begin
        air = NS.AirLayer{Float64}();
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.product_limited_rate!(ps, 20.0);
            @test true;

            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.product_limited_rate!(ps, air, 0.2);
            @test true;
        end;
    end;

    @testset "Photosynthesis colimitation" begin
        air = NS.AirLayer{Float64}();
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
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
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosystem_temperature_dependence!(ps, air, 298.15);
            PS.photosystem_electron_transport!(ps, 100.0, 20.0);
            PS.rubisco_limited_rate!(ps, 20.0);
            PS.light_limited_rate!(ps);
            PS.product_limited_rate!(ps, 20.0);
            PS.colimit_photosynthesis!(ps);
            PS.photosystem_coefficients!(ps, 100.0);
            @test true;
        end;

        # try out the qL based fluorescence model (CXVJP model only)
        for i in 1:3
            ps.trait = psts[i];
            ps.state = psss[i];
            ps.trait.FLM = NS.QLFluoscenceModel{Float64}();
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
        ps = NS.LeafPhotosystem{Float64}();
        psts = [NS.C3VJPTrait{Float64}(), NS.C4CLMTrait{Float64}(), NS.C4VJPTrait{Float64}(), NS.C3CytoTrait{Float64}()];
        psss = [NS.C3State{Float64}(), NS.C4State{Float64}(), NS.C4State{Float64}(), NS.C3State{Float64}()];
        for i in 1:4
            ps.trait = psts[i];
            ps.state = psss[i];
            PS.photosynthesis_only!(ps, air, 0.2, 100.0, 298.15);
            @test true;
        end;
    end;

    @testset "Derivatives (for stomatal models)" begin
        config = NS.SPACConfiguration(Float64);
        leaf = NS.Leaf(config);
        leaf.xylem.trait.k_max = 0.05;
        @test PS.∂R∂T(leaf) > 0;

        for td in [NS.Arrhenius{Float64}(T_REF = 298.15, VAL_REF = NaN , ΔHA = 46390.0),
                   NS.ArrheniusPeak{Float64}(T_REF = 298.15, VAL_REF = NaN , ΔHA = 46390.0, ΔHD = 150650.0, ΔSV = 490.0),
                   NS.Q10{Float64}(Q_10 = 1.4, T_REF = 298.15, VAL_REF = 0.0140/8760),
                   NS.Q10Peak{Float64}(Q_10 = 1.4, T_REF = 298.15, VAL_REF = 0.0140/8760, ΔHD = 150650.0, ΔSV = 490.0)]
            @test PS.∂R∂T(td, 1.0, 298.15) > 0;
        end;
    end;

    @testset "Leaf photosynthesis" begin
        config = NS.SPACConfiguration(Float64);
        leaf = NS.Leaf(config);
        leaf.xylem.trait.k_max = 0.05;
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar_sunlit .= 100.0;
        leaf.flux.auxil.ppar_shaded = 100.0;

        # respiration only
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = true);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(), 1.0; rd_only = true);
        @test true;

        # core model
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(), 1.0; rd_only = false);
        @test true;

        # optimality stomatal model
        leaf.flux.trait.stomatal_model = NS.WangSM{Float64}();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(); rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(); rd_only = false);
        @test true;

        # empirical stomatal model (beta on G1)
        leaf.flux.trait.stomatal_model = NS.BallBerrySM{Float64}();
        leaf.flux.trait.stomatal_model.β.PARAM_Y = NS.BetaParameterG1();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(); rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(); rd_only = false);
        @test true;

        # empirical stomatal model (beta on Vcmax)
        leaf.flux.trait.stomatal_model = NS.BallBerrySM{Float64}();
        leaf.flux.trait.stomatal_model.β.PARAM_Y = NS.BetaParameterVcmax();
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(); rd_only = false);
        @test true;
        PS.leaf_photosynthesis!(leaf, air, NS.PCO₂Mode(); rd_only = false);
        @test true;
    end;

    @testset "Plant photosynthesis" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);

        PS.plant_photosynthesis!(spac, NS.GCO₂Mode());
        @test true;
        PS.plant_photosynthesis!(spac, NS.PCO₂Mode());
        @test true;

        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        spac.canopy.structure.trait.lai = 0.0;
        PS.plant_photosynthesis!(spac, NS.GCO₂Mode());
        @test true;
        PS.plant_photosynthesis!(spac, NS.PCO₂Mode());
        @test true;

        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        spac.canopy.sun_geometry.state.sza = 90;
        PS.plant_photosynthesis!(spac, NS.GCO₂Mode());
        @test true;
        PS.plant_photosynthesis!(spac, NS.PCO₂Mode());
        @test true;
    end;
end;

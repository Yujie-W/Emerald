using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS
import Emerald.EmeraldLand.StomatalModels as SM
import Emerald.EmeraldLand.SPAC


@testset verbose = true "StomatalModels.jl" begin
    @testset "Empirical equations" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar_sunlit .= 100.0;
        leaf.flux.auxil.ppar_shaded = 100.0;
        leaf.flux.auxil.g_CO₂_shaded = 0.02;
        leaf.flux.auxil.g_CO₂_sunlit .= 0.02;
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = false);

        for sm in [NS.BallBerrySM{Float64}(), NS.GentineSM{Float64}(), NS.LeuningSM{Float64}(), NS.MedlynSM{Float64}()]
            leaf.flux.state.stomatal_model = sm;
            gsh = SM.empirical_equation(sm, leaf, air);
            gsl = SM.empirical_equation(sm, leaf, air, 1);
            @test gsh > 0.0;
            @test gsl > 0.0;
        end;
    end;

    @testset "Beta function" begin
        # compute the beta factor
        f(x) = x;
        @test SM.β_factor(f, NS.WeibullVC{Float64}(2, 5), 0.0) == 1.0;
        @test SM.β_factor(f, NS.WeibullVC{Float64}(2, 5), -1.0) < 1.0;
        @test SM.β_factor(f, NS.VanGenuchten{Float64}("Loam"), 0.0) == 1.0;
        @test SM.β_factor(f, NS.VanGenuchten{Float64}("Clay"), -1.0) < 1.0;
        @test SM.β_factor(f, 1.0) == 1.0;
        @test SM.β_factor(f, 0.5) < 1.0;

        # read the beta from stomatal models
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        SM.read_β(leaf);
        @test true;

        # set the beta factor based on stomatal models
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);

        # the function does not for optimality models
        SM.β_factor!(spac);
        for leaf in spac.LEAVES
            @test isnan(SM.read_β(leaf));
        end;

        # the function will set up the beta factor for empirical models
        for root in spac.ROOTS
            root.xylem.auxil.flow = 1.0;
        end;

        # BetaParameterKleaf
        for leaf in spac.LEAVES
            leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
            leaf.flux.state.stomatal_model.β.PARAM_X = NS.BetaParameterKleaf();
        end;
        SM.β_factor!(spac);
        for leaf in spac.LEAVES
            @test 0 < SM.read_β(leaf) <= 1;
        end;

        # BetaParameterKsoil
        for leaf in spac.LEAVES
            leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
            leaf.flux.state.stomatal_model.β.PARAM_X = NS.BetaParameterKsoil();
        end;
        SM.β_factor!(spac);
        for leaf in spac.LEAVES
            @test 0 < SM.read_β(leaf) <= 1;
        end;

        # BetaParameterPleaf
        for leaf in spac.LEAVES
            leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
            leaf.flux.state.stomatal_model.β.PARAM_X = NS.BetaParameterPleaf();
        end;
        SM.β_factor!(spac);
        for leaf in spac.LEAVES
            @test 0 < SM.read_β(leaf) <= 1;
        end;

        # BetaParameterPsoil
        for leaf in spac.LEAVES
            leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
            leaf.flux.state.stomatal_model.β.PARAM_X = NS.BetaParameterPsoil();
        end;
        SM.β_factor!(spac);
        for leaf in spac.LEAVES
            @test 0 < SM.read_β(leaf) <= 1;
        end;

        # BetaParameterΘ
        for leaf in spac.LEAVES
            leaf.flux.state.stomatal_model = NS.BallBerrySM{Float64}();
            leaf.flux.state.stomatal_model.β.PARAM_X = NS.BetaParameterΘ();
        end;
        SM.β_factor!(spac);
        for leaf in spac.LEAVES
            @test 0 < SM.read_β(leaf) <= 1;
        end;
    end;

    @testset "∂A∂E" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar_sunlit .= 100.0;
        leaf.flux.auxil.ppar_shaded = 100.0;
        leaf.flux.state.g_H₂O_s_shaded = 0.02;
        leaf.flux.state.g_H₂O_s_sunlit .= 0.02;
        SM.stomatal_conductance_profile!(leaf);
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = false);

        @test SM.∂A∂E(leaf, air) > 0;
        @test SM.∂A∂E(leaf, air, 1) > 0;
    end;

end;

using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS
import Emerald.EmeraldLand.PlantHydraulics as PH
import Emerald.EmeraldLand.StomatalModels as SM
import Emerald.EmeraldLand.SPAC


@testset verbose = true "StomatalModels.jl" begin
    config = NS.SPACConfiguration(Float64);
    spac = NS.BulkSPAC(config);

    @testset "Empirical equations" begin
        config = NS.SPACConfiguration(Float64);
        leaf = NS.CanopyLayer(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar .= 100.0;
        leaf.flux.auxil.g_CO₂ .= 0.02;
        PS.leaf_photosynthesis!(config, spac.cache, leaf, air, 1.0; rd_only = false);

        for sm in [NS.BallBerrySM{Float64}(), NS.GentineSM{Float64}(), NS.LeuningSM{Float64}(), NS.MedlynSM{Float64}()]
            leaf.flux.trait.stomatal_model = sm;
            gs = SM.empirical_equation(sm, leaf, air);
            @test all(gs .>= 0);
        end;
    end;

    @testset "Beta function" begin
        # compute the beta factor
        f(x) = x;
        @test SM.β_factor(f, 1.0) == 1.0;
        @test SM.β_factor(f, 0.5) < 1.0;

        # read the beta from stomatal models
        config = NS.SPACConfiguration(Float64);
        leaf = NS.CanopyLayer(config);
        SM.read_β(leaf);
        @test true;

        # set the beta factor based on stomatal models
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);

        # the function does not for optimality models
        SM.β_factor!(spac);
        for leaf in spac.plant.leaves
            @test isnan(SM.read_β(leaf));
        end;

        # the function will set up the beta factor for empirical models
        for root in spac.plant.roots
            root.xylem.auxil.flow = 1.0;
        end;

        # BetaParameterKleaf
        for param_x in [NS.BetaParameterKleaf(), NS.BetaParameterKsoil(), NS.BetaParameterPleaf(), NS.BetaParameterPsoil(), NS.BetaParameterΘ()]
            for leaf in spac.plant.leaves
                leaf.flux.trait.stomatal_model = NS.BallBerrySM{Float64}();
                leaf.flux.trait.stomatal_model.β.PARAM_X = param_x;
            end;
            SM.β_factor!(spac);
            for leaf in spac.plant.leaves
                @test 0 < SM.read_β(leaf) <= 1;
            end;
        end;
    end;

    @testset "∂A∂E" begin
        config = NS.SPACConfiguration(Float64);
        leaf = NS.CanopyLayer(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar .= 100.0;
        leaf.flux.state.g_H₂O_s .= 0.02;
        SPAC.substep_aux!(leaf);
        PS.leaf_photosynthesis!(config, spac.cache, leaf, air, 1.0; rd_only = false);
        SM.∂A∂E!(spac.cache, leaf, air);

        @test all(leaf.flux.auxil.∂A∂E .> 0);
    end;

    @testset "∂Θ∂E" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        SPAC.soil_plant_air_continuum!(config, spac, 10);
        leaf = spac.plant.leaves[1];
        air = spac.airs[1];

        for sm in [NS.WangSM{Float64}(), NS.AndereggSM{Float64}(), NS.EllerSM{Float64}(), NS.SperrySM{Float64}(), NS.Wang2SM{Float64}()]
            SM.∂Θ∂E!(spac.cache, sm, leaf, air);
            @test all(leaf.flux.auxil.∂Θ∂E .> 0);
        end;
    end;

    #=
    @testset "Nighttime model" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.state.g_H₂O_s_shaded = 0.02;
        leaf.flux.state.g_H₂O_s_sunlit .= 0.02;
        leaf.flux.auxil.ppar_sunlit .= 0;
        leaf.flux.auxil.ppar_shaded = 0;
        SPAC.substep_aux!(leaf);
        PS.leaf_photosynthesis!(leaf, air, 1.0; rd_only = false);
        PH.leaf_pressure_profile!(config, leaf, spac.cache, 0.0);

        @test SM.∂R∂E(leaf, air, 1.0) > 0;
        @test SM.∂Θₙ∂E(leaf, air) > 0;
    end;
    =#

    @testset "∂g∂t & ∂gₙ∂t" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        SPAC.soil_plant_air_continuum!(config, spac, 10);
        leaf = spac.plant.leaves[1];
        air = spac.airs[1];

        for sm in [NS.AndereggSM{Float64}(), NS.EllerSM{Float64}(), NS.SperrySM{Float64}(), NS.WangSM{Float64}(), NS.Wang2SM{Float64}()]
            SM.∂g∂t!(spac.cache, sm, leaf, air);
            @test all(leaf.flux.auxil.∂g∂t .> 0);
        end;

        for sm in [NS.BallBerrySM{Float64}(), NS.GentineSM{Float64}(), NS.LeuningSM{Float64}(), NS.MedlynSM{Float64}()]
            leaf.flux.auxil.β = 0.9;
            sm.β.PARAM_Y = NS.BetaParameterG1();
            SM.∂g∂t!(spac.cache, sm, leaf, air);
            @test all(leaf.flux.auxil.∂g∂t .> 0);

            sm.β.PARAM_Y = NS.BetaParameterVcmax();
            SM.∂g∂t!(spac.cache, sm, leaf, air);
            @test all(leaf.flux.auxil.∂g∂t .> 0);
        end;

        # ∂gₙ∂t is only valid for WangSM
        leaf.flux.state.g_H₂O_s .= 0.1;
        SPAC.substep_aux!(leaf);
        @test SM.∂gₙ∂t(leaf, air, 1.0) < 0;
    end;

    @testset "Stomatal limits" begin
        config = NS.SPACConfiguration(Float64);
        leaf = NS.CanopyLayer(config);

        leaf.flux.state.g_H₂O_s .= 0;
        SM.limit_stomatal_conductance!(leaf);
        @test all(leaf.flux.state.g_H₂O_s .> 0);

        leaf.flux.state.g_H₂O_s .= 1;
        SM.limit_stomatal_conductance!(leaf);
        @test all(leaf.flux.state.g_H₂O_s .< 1);
    end;

    @testset "Stomatal profiles" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        for leaf in spac.plant.leaves
            leaf.flux.auxil.ppar .= 200.0;
            leaf.flux.state.g_H₂O_s .= 0;
            SM.limit_stomatal_conductance!(leaf);
        end;
        PS.plant_photosynthesis!(config, spac);
        SM.stomatal_conductance_profile!(spac);

        for leaf in spac.plant.leaves
            @test all(leaf.flux.auxil.∂g∂t .> 0);
        end;
    end;

    @testset "Stomatal budgets" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        gss = [];
        for leaf in spac.plant.leaves
            leaf.flux.auxil.ppar .= 200.0;
            leaf.flux.state.g_H₂O_s .= 0;
            SM.limit_stomatal_conductance!(leaf);
            push!(gss, deepcopy(leaf.flux.state.g_H₂O_s));
        end;
        PS.plant_photosynthesis!(config, spac);
        SM.stomatal_conductance_profile!(spac);
        SM.stomatal_conductance!(spac, 1.0);

        for i in eachindex(spac.plant.leaves)
            leaf = spac.plant.leaves[i];
            @test all(leaf.flux.auxil.∂g∂t .> 0);

            @test all(leaf.flux.state.g_H₂O_s .== gss[i] .+ leaf.flux.auxil.∂g∂t .* 1.0);
        end;
    end;
end;

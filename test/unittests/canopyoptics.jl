using Test
import Emerald.EmeraldLand.CanopyOptics as CO
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.SPAC


@testset verbose = true "Canopy Optics Model" begin
    @testset "Canopy extinction coefficient" begin
        sza = 40.0;
        lias = collect(Float64, 0:5:90);

        ks = CO.extinction_coefficient.(sza, lias);
        @test all(ks .>= 0);

        kd = CO.extinction_coefficient.(lias);
        @test all(kd .>= 0);
    end;

    @testset "Soil albedo" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);

        spac.SOIL_BULK.auxil._θ = -1;
        config.α_CLM = true;
        config.α_FITTING = false;
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.SOIL_BULK.auxil.ρ_sw .< 1);

        spac.SOIL_BULK.auxil._θ = -1;
        config.α_CLM = false;
        config.α_FITTING = false;
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.SOIL_BULK.auxil.ρ_sw .< 1);

        spac.SOIL_BULK.auxil._θ = -1;
        config.α_CLM = false;
        config.α_FITTING = true;
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.SOIL_BULK.auxil.ρ_sw .< 1);
    end;

    @testset "Canopy structure" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        CO.canopy_structure!(config, spac);

        @test 0 <= spac.CANOPY.structure.auxil.bf <= 1;
        @test spac.CANOPY.structure.auxil.ddb >= 0;
        @test spac.CANOPY.structure.auxil.ddf >= 0;
        @test all(0 .< spac.CANOPY.structure.auxil.ρ_lw_layer .< 1);
        @test all(0 .< spac.CANOPY.structure.auxil.τ_lw_layer .< 1);
        @test all(0 .< spac.CANOPY.structure.auxil.ϵ_lw_layer .< 1);
        @test all(0 .< spac.CANOPY.structure.auxil.ρ_lw .< 1);
        @test all(0 .< spac.CANOPY.structure.auxil.τ_lw .< 1);
    end;

    @testset "Canopy sun geometry" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry!(config, spac);

        @test spac.CANOPY.sun_geometry.auxil.ks >= 0;
        @test spac.CANOPY.sun_geometry.auxil.sdb >= 0;
        @test spac.CANOPY.sun_geometry.auxil.sdf >= 0;
        @test 0 < spac.CANOPY.structure.auxil.ci <= 1;
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.p_sunlit .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.ρ_sd_layer .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.ρ_dd_layer .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.τ_ss_layer .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.τ_sd_layer .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.τ_dd_layer .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.ρ_sd .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.ρ_dd .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.τ_sd .< 1);
        @test all(0 .< spac.CANOPY.sun_geometry.auxil.τ_dd .< 1);
    end;

    @testset "Canopy sensor geometry" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.sensor_geometry!(config, spac);

        @test spac.CANOPY.sensor_geometry.auxil.ko >= 0;
        @test spac.CANOPY.sensor_geometry.auxil.dob >= 0;
        @test spac.CANOPY.sensor_geometry.auxil.dof >= 0;
        @test spac.CANOPY.sensor_geometry.auxil.sob >= 0;
        @test spac.CANOPY.sensor_geometry.auxil.sof >= 0;
        @test all(0 .< spac.CANOPY.sensor_geometry.auxil.p_sensor .< 1);
        @test 0 < spac.CANOPY.sensor_geometry.auxil.p_sensor_soil < 1
        @test all(0 .< spac.CANOPY.sensor_geometry.auxil.p_sun_sensor .< 1);
    end;

    @testset "Shortwave radiation" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.shortwave_radiation!(config, spac);

        @test all(spac.CANOPY.sun_geometry.auxil.e_dirꜜ .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_difꜜ .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_difꜛ .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_net_dir .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_net_dif .> 0);
        @test all(spac.SOIL_BULK.auxil.e_net_dir .> 0);
        @test all(spac.SOIL_BULK.auxil.e_net_dif .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.r_net_sw .> 0);
        @test all(spac.SOIL_BULK.auxil.r_net_sw .> 0);
        for leaf in spac.LEAVES
            @test all(leaf.flux.auxil.apar_shaded .> 0);
            @test all(leaf.flux.auxil.apar_sunlit .> 0);
            @test all(leaf.flux.auxil.ppar_shaded .> 0);
            @test all(leaf.flux.auxil.ppar_sunlit .> 0);
        end;
    end;

    @testset "Longwave radiation" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.longwave_radiation!(config, spac);

        @test all(spac.CANOPY.structure.auxil.lw_layer .> 0);
        @test all(spac.CANOPY.structure.auxil.emitꜜ .> 0);
        @test all(spac.CANOPY.structure.auxil.emitꜛ .> 0);
        @test all(spac.CANOPY.structure.auxil.lwꜜ .> 0);
        @test all(spac.CANOPY.structure.auxil.lwꜛ .>= 0);
        @test all(!isnan, spac.CANOPY.structure.auxil.r_net_lw);
        @test all(!isnan, spac.SOIL_BULK.auxil.r_net_lw);
    end;

    @testset "Canopy reflection spectrum" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.shortwave_radiation!(config, spac);
        CO.sensor_geometry!(config, spac);
        CO.reflection_spectrum!(config, spac);

        @test all(spac.CANOPY.sensor_geometry.auxil.e_sensor_layer .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.e_sensor .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.reflectance .> 0);
    end;

    @testset "Canopy fluorescence spectrum" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        for leaf in spac.LEAVES
            leaf.flux.auxil.ϕ_f_shaded = 0.01;
            leaf.flux.auxil.ϕ_f_sunlit .= 0.01;
        end;
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.shortwave_radiation!(config, spac);
        CO.sensor_geometry!(config, spac);
        CO.fluorescence_spectrum!(config, spac);

        @test all(spac.CANOPY.sun_geometry.auxil.e_sif_chl .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_sunlit .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_shaded .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_scattered .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜜ_layer .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜛ_layer .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜜ_emit .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜛ_emit[:,1:end-1] .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜛ_emit[:,end] .== 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜜ[:,1] .== 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜜ[:,2:end] .> 0);
        @test all(spac.CANOPY.sun_geometry.auxil.e_sifꜛ .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_obs_sunlit .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_obs_shaded .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_obs_scattered .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_obs_soil .> 0);
        @test all(spac.CANOPY.sensor_geometry.auxil.sif_obs .> 0);
    end;

end;

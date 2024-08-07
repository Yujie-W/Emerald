using Test
import Emerald.EmeraldLand.CanopyOptics as CO
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.SPAC


@testset verbose = true "CanopyOptics.jl" begin
    @testset "Canopy extinction coefficient" begin
        sza = 40.0;
        lias = collect(Float64, 0:5:90);

        ks = CO.extinction_coefficient.(sza, lias);
        @test all(ks .>= 0);

        kd = CO.extinction_coefficient.(lias);
        @test all(kd .>= 0);
    end;

    @testset "Soil albedo" begin
        config = NS.SPACConfiguration(Float64);

        config.SOIL_ALBEDO = NS.SoilAlbedoHyperspectralCLM();
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.soil_bulk.auxil.ρ_sw .< 1);

        config.SOIL_ALBEDO = NS.SoilAlbedoBroadbandCLM();
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.soil_bulk.auxil.ρ_sw .< 1);

        config.SOIL_ALBEDO = NS.SoilAlbedoHyperspectralCLIMA();
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.soil_bulk.auxil.ρ_sw .< 1);

        config.SOIL_ALBEDO = NS.SoilAlbedoBroadbandCLIMA();
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        @test all(0 .< spac.soil_bulk.auxil.ρ_sw .< 1);
    end;

    @testset "Canopy structure" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);

        @test spac.canopy.structure.t_aux.ddb >= 0;
        @test spac.canopy.structure.t_aux.ddf >= 0;
        @test all(0 .< spac.canopy.structure.auxil.ρ_dd_layer .< 1);
        @test all(0 .< spac.canopy.structure.auxil.τ_dd_layer .< 1);
        @test all(0 .< spac.canopy.structure.auxil.ρ_dd .< 1);
        @test all(0 .< spac.canopy.structure.auxil.τ_dd .< 1);
        @test all(0 .< spac.canopy.structure.auxil.ρ_lw_layer .< 1);
        @test all(0 .< spac.canopy.structure.auxil.τ_lw_layer .< 1);
        @test all(0 .< spac.canopy.structure.auxil.ϵ_lw_layer .< 1);
        @test all(0 .< spac.canopy.structure.auxil.ρ_lw .< 1);
        @test all(0 .< spac.canopy.structure.auxil.τ_lw .< 1);
    end;

    @testset "Canopy sun geometry" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry_aux!(config, spac);
        CO.sun_geometry!(config, spac);

        @test spac.canopy.sun_geometry.s_aux.ks >= 0;
        @test spac.canopy.sun_geometry.s_aux.sdb >= 0;
        @test spac.canopy.sun_geometry.s_aux.sdf >= 0;
        @test 0 < spac.canopy.structure.trait.ci <= 1;
        @test all(0 .< spac.canopy.sun_geometry.s_aux.p_sunlit .< 1);
        @test all(0 .< spac.canopy.sun_geometry.auxil.ρ_sd_layer .< 1);
        @test all(0 .< spac.canopy.sun_geometry.auxil.τ_ss_layer .< 1);
        @test all(0 .< spac.canopy.sun_geometry.auxil.τ_sd_layer .< 1);
        @test all(0 .< spac.canopy.sun_geometry.auxil.ρ_sd .< 1);
        @test all(0 .< spac.canopy.sun_geometry.auxil.τ_sd .< 1);
    end;

    @testset "Canopy sensor geometry" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry_aux!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.sensor_geometry_aux!(config, spac);
        CO.sensor_geometry!(config, spac);

        @test spac.canopy.sensor_geometry.s_aux.ko >= 0;
        @test spac.canopy.sensor_geometry.s_aux.dob >= 0;
        @test spac.canopy.sensor_geometry.s_aux.dof >= 0;
        @test spac.canopy.sensor_geometry.s_aux.sob >= 0;
        @test spac.canopy.sensor_geometry.s_aux.sof >= 0;
        @test all(0 .< spac.canopy.sensor_geometry.s_aux.p_sensor .< 1);
        @test 0 < spac.canopy.sensor_geometry.s_aux.p_sensor_soil < 1;
        @test all(0 .< spac.canopy.sensor_geometry.s_aux.p_sun_sensor .< 1);
        @test all(spac.canopy.sensor_geometry.s_aux.p_sun_sensor .< spac.canopy.sun_geometry.s_aux.p_sunlit);
        @test all(spac.canopy.sensor_geometry.s_aux.p_sun_sensor .< spac.canopy.sensor_geometry.s_aux.p_sensor);

        config = NS.SPACConfiguration(Float64);
        tpac = NS.BulkSPAC(config);
        SPAC.prescribe_traits!(config, tpac; ci = 0.5);
        SPAC.initialize_spac!(config, tpac);
        CO.soil_albedo!(config, tpac);
        CO.canopy_structure!(config, tpac);
        CO.sun_geometry_aux!(config, tpac);
        CO.sun_geometry!(config, tpac);
        CO.sensor_geometry_aux!(config, tpac);
        CO.sensor_geometry!(config, tpac);
        @test all(0 .< tpac.canopy.sun_geometry.s_aux.p_sunlit .< 1);
        @test all(tpac.canopy.sensor_geometry.s_aux.p_sun_sensor .< tpac.canopy.sun_geometry.s_aux.p_sunlit);
        @test all(tpac.canopy.sensor_geometry.s_aux.p_sun_sensor .< tpac.canopy.sensor_geometry.s_aux.p_sensor);
    end;

    @testset "Shortwave radiation" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry_aux!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.shortwave_radiation!(config, spac);

        @test all(spac.canopy.sun_geometry.auxil.e_dirꜜ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜜ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜛ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dir .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dif .> 0);
        @test all(spac.soil_bulk.auxil.e_net_dir .> 0);
        @test all(spac.soil_bulk.auxil.e_net_dif .> 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_leaf .> 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_stem .> 0);
        @test all(spac.soil_bulk.auxil.r_net_sw .> 0);
        for leaf in spac.plant.leaves
            @test all(leaf.flux.auxil.ppar .> 0);
        end;
    end;

    @testset "Longwave radiation" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.longwave_radiation!(spac);

        @test all(spac.canopy.structure.auxil.lw_layer .> 0);
        @test all(spac.canopy.structure.auxil.emitꜜ .> 0);
        @test all(spac.canopy.structure.auxil.emitꜛ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜜ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜛ .>= 0);
        @test all(!isnan, spac.canopy.structure.auxil.r_net_lw_leaf);
        @test all(!isnan, spac.canopy.structure.auxil.r_net_lw_stem);
        @test all(!isnan, spac.soil_bulk.auxil.r_net_lw);
    end;

    @testset "Canopy reflection spectrum" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry_aux!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.shortwave_radiation!(config, spac);
        CO.sensor_geometry_aux!(config, spac);
        CO.sensor_geometry!(config, spac);
        CO.reflection_spectrum!(config, spac);

        @test all(spac.canopy.sensor_geometry.auxil.e_sensor_layer .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.e_sensor .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.reflectance .> 0);
    end;

    @testset "Canopy fluorescence spectrum" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);
        for leaf in spac.plant.leaves
            leaf.photosystem.auxil.ϕ_f .= 0.01;
            leaf.photosystem.auxil.ϕ_f1 .= 0.01;
            leaf.photosystem.auxil.ϕ_f2 .= 0.01;
        end;
        CO.soil_albedo!(config, spac);
        CO.canopy_structure!(config, spac);
        CO.sun_geometry_aux!(config, spac);
        CO.sun_geometry!(config, spac);
        CO.shortwave_radiation!(config, spac);
        CO.sensor_geometry_aux!(config, spac);
        CO.sensor_geometry!(config, spac);
        CO.fluorescence_spectrum!(config, spac);

        @test all(spac.canopy.sun_geometry.auxil.e_sif_chl .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_sunlit .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_shaded .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_scattered .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜜ_layer .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜛ_layer .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜜ_emit .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜛ_emit[:,1:end-1] .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜛ_emit[:,end] .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜜ[:,1] .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜜ[:,2:end] .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_sifꜛ .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_obs_sunlit .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_obs_shaded .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_obs_scattered .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_obs_soil .> 0);
        @test all(spac.canopy.sensor_geometry.auxil.sif_obs .> 0);
    end;

    @testset "Canopy radiation" begin
        config = NS.SPACConfiguration(Float64);
        spac = NS.BulkSPAC(config);
        SPAC.initialize_spac!(config, spac);

        # SZA < 90
        spac.canopy.sun_geometry.state.sza = 30;
        CO.canopy_radiation!(config, spac);
        @test all(spac.canopy.sun_geometry.auxil.e_dirꜜ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜜ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜛ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dir .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dif .> 0);
        @test all(spac.soil_bulk.auxil.e_net_dir .> 0);
        @test all(spac.soil_bulk.auxil.e_net_dif .> 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_leaf .> 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_stem .> 0);
        @test all(spac.soil_bulk.auxil.r_net_sw .> 0);
        for leaf in spac.plant.leaves
            @test all(leaf.flux.auxil.ppar .> 0);
        end;
        @test all(spac.canopy.structure.auxil.lw_layer .> 0);
        @test all(spac.canopy.structure.auxil.emitꜜ .> 0);
        @test all(spac.canopy.structure.auxil.emitꜛ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜜ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜛ .> 0);
        @test all(!isnan, spac.canopy.structure.auxil.r_net_lw_leaf);
        @test all(!isnan, spac.canopy.structure.auxil.r_net_lw_stem);
        @test all(!isnan, spac.soil_bulk.auxil.r_net_lw);

        # SZA > = 90
        spac.canopy.sun_geometry.state.sza = 90;
        CO.canopy_radiation!(config, spac);
        @test all(spac.canopy.sun_geometry.auxil.e_dirꜜ .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜜ .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜛ .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dir .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dif .== 0);
        @test all(spac.soil_bulk.auxil.e_net_dir .== 0);
        @test all(spac.soil_bulk.auxil.e_net_dif .== 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_leaf .== 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_stem .== 0);
        @test all(spac.soil_bulk.auxil.r_net_sw .== 0);
        for leaf in spac.plant.leaves
            @test all(leaf.flux.auxil.ppar .== 0);
        end;
        @test all(spac.canopy.structure.auxil.lw_layer .> 0);
        @test all(spac.canopy.structure.auxil.emitꜜ .> 0);
        @test all(spac.canopy.structure.auxil.emitꜛ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜜ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜛ .>= 0);
        @test all(!isnan, spac.canopy.structure.auxil.r_net_lw_leaf);
        @test all(!isnan, spac.canopy.structure.auxil.r_net_lw_stem);
        @test all(!isnan, spac.soil_bulk.auxil.r_net_lw);

        # LAI <= 0 and SAI <= 0
        spac.canopy.sun_geometry.state.sza = 30;
        spac.canopy.structure.trait.lai = 0;
        spac.canopy.structure.trait.sai = 0;
        CO.canopy_radiation!(config, spac);
        @test all(spac.canopy.sun_geometry.auxil.e_dirꜜ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜜ .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜛ[:,1:end-1] .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_difꜛ[:,end] .> 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dir .== 0);
        @test all(spac.canopy.sun_geometry.auxil.e_net_dif .== 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_leaf .== 0);
        @test all(spac.canopy.sun_geometry.auxil.r_net_sw_stem .== 0);
        @test all(spac.soil_bulk.auxil.e_net_dir .> 0);
        @test all(spac.soil_bulk.auxil.e_net_dif .> 0);
        @test all(spac.soil_bulk.auxil.r_net_sw .> 0);
        for leaf in spac.plant.leaves
            @test all(leaf.flux.auxil.ppar .== 0);
        end;
        @test all(spac.canopy.structure.auxil.lw_layer .== 0);
        @test all(spac.canopy.structure.auxil.emitꜜ .== 0);
        @test all(spac.canopy.structure.auxil.emitꜛ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜜ .> 0);
        @test all(spac.canopy.structure.auxil.lwꜛ .> 0);
        @test all(spac.canopy.structure.auxil.r_net_lw_leaf .== 0);
        @test all(spac.canopy.structure.auxil.r_net_lw_stem .== 0);
        @test all(!isnan, spac.soil_bulk.auxil.r_net_lw);
    end;
end;

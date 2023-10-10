using Test
import Emerald.EmeraldLand.CanopyOptics as CO
import Emerald.EmeraldLand.Namespace as NS


@testset verbose = true "Canopy Optics Model" begin
    @testset "Canopy extinction coefficient" begin
        sza = 40.0;
        lias = collect(Float64, 0:5:90);

        ks = CO.extinction_coefficient.(sza, lias);
        @test all(ks .>= 0);

        kd = CO.extinction_coefficient.(lias);
        @test all(kd .>= 0);
    end;

    @testset "Canopy structure" begin
        config = NS.SPACConfiguration{Float64}();
        can = NS.MultiLayerCanopy(config);
        CO.canopy_structure!(config, can);
        @test 0 <= can.structure.auxil.bf <= 1;
    end;

    @testset "Canopy sun geometry" begin
        config = NS.SPACConfiguration{Float64}();
        can = NS.MultiLayerCanopy(config);
        CO.canopy_structure!(config, can);
        CO.sun_geometry!(config, can);

        @test can.sun_geometry.auxil.ks >= 0;
        @test can.sun_geometry.auxil.ddb >= 0;
        @test can.sun_geometry.auxil.ddf >= 0;
        @test can.sun_geometry.auxil.sdb >= 0;
        @test can.sun_geometry.auxil.sdf >= 0;
        @test 0 < can.structure.auxil.ci <= 1;
        @test all(0 .< can.sun_geometry.auxil.ps .<= 1);
        @test all(0 .< can.sun_geometry.auxil.p_sunlit .< 1);
    end;

    @testset "Canopy sensor geometry" begin
        config = NS.SPACConfiguration{Float64}();
        can = NS.MultiLayerCanopy(config);
        CO.canopy_structure!(config, can);
        CO.sun_geometry!(config, can);
        CO.sensor_geometry!(config, can);

        @test can.sensor_geometry.auxil.ko >= 0;
        @test can.sensor_geometry.auxil.dob >= 0;
        @test can.sensor_geometry.auxil.dof >= 0;
        @test can.sensor_geometry.auxil.sob >= 0;
        @test can.sensor_geometry.auxil.sof >= 0;
        @test all(0 .< can.sensor_geometry.auxil.po .<= 1);
        @test all(0 .< can.sensor_geometry.auxil.pso .< 1);
    end;

    @testset "Soil albedo" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);

        spac.SOIL_BULK.auxil._θ = -1;
        config.α_CLM = true;
        config.α_FITTING = false;
        CO.soil_albedo!(config, spac.SOIL_BULK, spac.SOILS[1]);
        @test all(0 .< spac.SOIL_BULK.auxil.ρ_sw .< 1);

        spac.SOIL_BULK.auxil._θ = -1;
        config.α_CLM = false;
        config.α_FITTING = false;
        CO.soil_albedo!(config, spac.SOIL_BULK, spac.SOILS[1]);
        @test all(0 .< spac.SOIL_BULK.auxil.ρ_sw .< 1);

        spac.SOIL_BULK.auxil._θ = -1;
        config.α_CLM = false;
        config.α_FITTING = true;
        CO.soil_albedo!(config, spac.SOIL_BULK, spac.SOILS[1]);
        @test all(0 .< spac.SOIL_BULK.auxil.ρ_sw .< 1);
    end;

end;

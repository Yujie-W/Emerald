import Emerald.Land as ELAND
using DataFrames
using Test


@testset "Emerald Land" verbose = true begin
    @testset "Testing Mode with FvCB model" begin
        settings = ELAND.land_model_settings(mode = "testing");
        settings["C3_MODEL"] = "FvCB";
        @test true;
        df = ELAND.simulation!(settings, 31.86389, 117.28083, 2019);
        @test typeof(df) == DataFrame;
    end;

    @testset "Testing Mode with J3B model" begin
        settings = ELAND.land_model_settings(mode = "testing");
        settings["C3_MODEL"] = "J3B";
        @test true;
        df = ELAND.simulation!(settings, 31.86389, 117.28083, 2019);
        @test typeof(df) == DataFrame;
    end;
end;

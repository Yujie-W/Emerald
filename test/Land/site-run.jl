using DataFrames
using Test

import Emerald.Land as ELAND


@testset "Emerald Land" verbose = true begin
    @testset "Testing Mode with C3-Jmax-KN model" begin
        settings = ELAND.land_model_settings(mode = "testing");
        settings["C3_MODEL"] = "Jmax";
        settings["C3_ΦF_MODEL"] = "KN";
        @test true;
        nt = ELAND.simulation!(settings, 31.86389, 117.28083, 2019);
        @test typeof(nt) <: NamedTuple;
    end;

    @testset "Testing Mode with C3-Jmax-QL model" begin
        settings = ELAND.land_model_settings(mode = "testing");
        settings["C3_MODEL"] = "Jmax";
        settings["C3_ΦF_MODEL"] = "QL";
        @test true;
        nt = ELAND.simulation!(settings, 31.86389, 117.28083, 2019);
        @test typeof(nt) <: NamedTuple;
    end;

    @testset "Testing Mode with C3-Vqmax model" begin
        settings = ELAND.land_model_settings(mode = "testing");
        settings["C3_MODEL"] = "Vqmax";
        @test true;
        nt = ELAND.simulation!(settings, 31.86389, 117.28083, 2019);
        @test typeof(nt) <: NamedTuple;
    end;

    @testset "Testing Mode with C4-Vcmax-KN model" begin
        settings = ELAND.land_model_settings(mode = "testing");
        settings["C4_MODEL"] = "Vcmax";
        settings["C4_ΦF_MODEL"] = "KN";
        @test true;
        nt = ELAND.simulation!(settings, 31.86389, 117.28083, 2019; c3c4 = "C4");
        @test typeof(nt) <: NamedTuple;
    end;
end;

using PkgUtility.DataIO: read_csv
using Test

import Emerald.Namespace as ENS
import Emerald.ResearchTools as ERT


@testset "Emerald ResearchTools" verbose = true begin
    df3 = read_csv(joinpath(@__DIR__, "../..", "data/examples", "C3-ACi.csv"));
    df4 = read_csv(joinpath(@__DIR__, "../..", "data/examples", "C4-ACi.csv"));
    df3.T_LEAF .+= 273.15;  # convert to Kelvin
    df4.T_LEAF .+= 273.15;  # convert to Kelvin

    @testset "C3 Jmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C3_AC_METHOD = ENS.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = ENS.AjMethodC3JmaxPi();
        config.METHODS.C3_AP_METHOD = ENS.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = ENS.ColimitJCLM(Float64);
        config.METHODS.FLUORESCENCE_METHOD_C3 = ENS.KNFluorescenceModel{Float64}();
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25", "Rd25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25", "Γstar25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;

    @testset "C3 Vqmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C3_AC_METHOD = ENS.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = ENS.AjMethodC3VqmaxPi();
        config.METHODS.C3_AP_METHOD = ENS.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = ENS.SerialColimit();
        config.METHODS.FLUORESCENCE_METHOD_C3 = ENS.CytochromeFluorescenceModel();
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "b₆f"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "b₆f", "Rd25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "b₆f", "Γstar25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;

    @testset "C4 Vcmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C4_AP_METHOD = ENS.ApMethodC4VcmaxPi();
        config.METHODS.FLUORESCENCE_METHOD_C4 = ENS.KNFluorescenceModel{Float64}();
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;

    @testset "C4 Vpmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C4_AP_METHOD = ENS.ApMethodC4VpmaxPi();
        config.METHODS.FLUORESCENCE_METHOD_C4 = ENS.KNFluorescenceModel{Float64}();
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25", "Vpmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25", "Vpmax25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;
end;

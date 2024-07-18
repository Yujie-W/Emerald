# This is an tutorial about how to turn off the reflectance and SIF features to speed up the simulation
using Emerald
using Test

@testset "Turn off Reflectance and/or SIF" begin
    FT = Float64;

    # default configuration
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @time EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @test true;

    # turn off reflectance and SIF
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    config.ENABLE_REF = false;
    config.ENABLE_SIF = false;
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @time EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @test true;
end;

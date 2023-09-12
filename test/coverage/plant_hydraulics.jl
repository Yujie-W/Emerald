@testset verbose = true "PlantHydraulics" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.MultiLayerSPAC(config);
    EmeraldLand.SPAC.initialize!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    EmeraldLand.SPAC.update!(config, spac; swcs = (0.08,0.09,0.10,0.2,0.3));
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    EmeraldLand.SPAC.update!(config, spac; swcs = (0.08,0.08,0.08,0.2,0.3));
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    EmeraldLand.SPAC.update!(config, spac; swcs = (0.08,0.08,0.08,0.08,0.08));
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    @test true;
end

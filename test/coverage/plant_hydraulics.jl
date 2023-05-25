@testset verbose = true "PlantHydraulics" begin
    FT = Float64;
    spac = EmeraldLand.Namespace.MultiLayerSPAC{FT}();
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    EmeraldLand.SPAC.initialize!(spac, config);
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    EmeraldLand.SPAC.update!(spac, config; swcs = (0.08,0.09,0.10,0.2,0.3));
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    EmeraldLand.SPAC.update!(spac, config; swcs = (0.08,0.08,0.08,0.2,0.3));
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    EmeraldLand.SPAC.update!(spac, config; swcs = (0.08,0.08,0.08,0.08,0.08));
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    @test true;
end

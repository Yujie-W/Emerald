@testset verbose = true "PlantHydraulics" begin
    FT = Float64;
    config = EmeraldCore.Namespace.SPACConfiguration{FT}(20);
    spac = EmeraldCore.Namespace.MultiLayerSPAC(config);
    EmeraldCore.SPAC.initialize!(spac, config);
    EmeraldCore.SPAC.spac!(spac, config, FT(1));

    EmeraldCore.SPAC.update!(spac, config; swcs = (0.08,0.09,0.10,0.2,0.3));
    EmeraldCore.SPAC.spac!(spac, config, FT(1));

    EmeraldCore.SPAC.update!(spac, config; swcs = (0.08,0.08,0.08,0.2,0.3));
    EmeraldCore.SPAC.spac!(spac, config, FT(1));

    EmeraldCore.SPAC.update!(spac, config; swcs = (0.08,0.08,0.08,0.08,0.08));
    EmeraldCore.SPAC.spac!(spac, config, FT(1));

    @test true;
end

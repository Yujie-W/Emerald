@testset verbose = true "EmeraldCore" begin
    FT = Float64;
    gcf = EmeraldCore.Namespace.GeneralConfiguration();
    dims = EmeraldCore.Namespace.SPACDimension(gcf);
    config = EmeraldCore.Namespace.SPACConfiguration{FT,dims}(gcf);
    spac = EmeraldCore.Namespace.MultiLayerSPAC(config);
    EmeraldCore.SPAC.initialize!(spac, config);
    EmeraldCore.SPAC.spac!(spac, config, FT(1));
    @test true;
end

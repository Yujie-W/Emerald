@testset verbose = true "EmeraldCore" begin
    FT = Float64;
    config = EmeraldCore.Namespace.SPACConfiguration{FT}();
    spac = EmeraldCore.Namespace.MultiLayerSPAC{FT}();
    EmeraldCore.SPAC.initialize!(spac);
    EmeraldCore.SPAC.spac!(spac, config, FT(1));
    @test true;
end

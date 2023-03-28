@testset verbose = true "PlantHydraulics" begin
    FT = Float64;
    spac = EmeraldCore.Namespace.MultiLayerSPAC{FT}();
    EmeraldCore.SPAC.initialize!(spac);
    @time EmeraldCore.SPAC.spac!(spac, FT(1));

    EmeraldCore.SPAC.update!(spac; swcs = (0.08,0.09,0.10,0.2,0.3));
    @time EmeraldCore.SPAC.spac!(spac, FT(1));

    EmeraldCore.SPAC.update!(spac; swcs = (0.08,0.08,0.08,0.2,0.3));
    @time EmeraldCore.SPAC.spac!(spac, FT(1));

    EmeraldCore.SPAC.update!(spac; swcs = (0.08,0.08,0.08,0.08,0.08));
    @time EmeraldCore.SPAC.spac!(spac, FT(1));

    @test true;
end

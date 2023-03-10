@testset verbose = true "EmeraldCore" begin
    FT = Float64;
    spac = EmeraldCore.Namespace.MonoMLTreeSPAC{FT}();
    EmeraldCore.SPAC.initialize!(spac);
    EmeraldCore.SPAC.spac!(spac, FT(1));
    @test true;
end

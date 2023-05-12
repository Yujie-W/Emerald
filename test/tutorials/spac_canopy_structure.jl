@testset "Modify Canopy Structural Parameters" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.MultiLayerSPAC{FT}();
    EmeraldLand.SPAC.initialize!(spac, config);
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    # Changing canopy structure may result in changes in other parameters, such as Vcmax profile.
    # Thus, it is not recommended to modify those parametes manually unless otherwise told to by developers.
    # It is highly recommended to use our embedded function update! to modify canopy structure.
    # Currently, function update! supports the modification of leaf area index and clumping index.
    EmeraldLand.SPAC.update!(spac, config; lai = 3, ci = 0.8);
    @test true;
end

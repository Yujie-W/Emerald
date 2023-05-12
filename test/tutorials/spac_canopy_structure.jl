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

    # Leaf inclination angle distribution is stored as a vector P_INCL in field CANOPY.
    # By default, P_INCL is a uniform distribution from 0° to 90° per 10°.
    # To change P_INCL, we use the VerhoefLIDF model to update the distribution function.
    # For example, (A,B) =
    #     (0,0) gives uniform distribution,
    #     (-0.35,-0.15) gives spherical distribution,
    #     (-1,0) gives erectophile distribution,
    #     (1,0) gives planophile distribution,
    #     (0,-1) gives plagiophile distribution, and
    #     (0,1) gives extremophile distribution.
    spac.CANOPY.LIDF.A = -0.35;
    spac.CANOPY.LIDF.B = -0.15;
    EmeraldLand.CanopyOptics.inclination_angles!(spac.CANOPY, spac.CANOPY.LIDF);
    EmeraldLand.SPAC.spac!(spac, config, FT(1));
    @test true;
end

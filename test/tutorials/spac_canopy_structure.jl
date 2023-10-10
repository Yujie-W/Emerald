using Emerald;
using Test;


@testset "Modify Canopy Structural Parameters" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.MultiLayerSPAC(config);
    EmeraldLand.SPAC.initialize!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    # Changing canopy structure may result in changes in other parameters, such as Vcmax profile.
    # Thus, it is not recommended to modify those parametes manually unless otherwise told to by developers.
    # It is highly recommended to use our embedded function update! to modify canopy structure.
    # Currently, function update! supports the modification of leaf area index and clumping index.
    EmeraldLand.SPAC.update!(config, spac; lai = 3, ci = 0.8);
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
    spac.CANOPY.structure.state.LIDF.A = -0.35;
    spac.CANOPY.structure.state.LIDF.B = -0.15;
    EmeraldLand.CanopyOptics.inclination_angles!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;

    # We also support beta function distribution. For example (A,B) =
    #     (1.000,1.000) gives uniform distribution,
    #     (1.930,1.101) gives spherical distribution,
    #     (2.770,1.172) gives erectophile distribution,
    #     (1.172,2.770) gives planophile distribution,
    #     (3.326,3.326) gives plagiophile distribution, and
    #     (0.433,0.433) gives extremophile distribution.
    spac.CANOPY.structure.state.LIDF = EmeraldLand.Namespace.BetaLIDF{FT}();
    spac.CANOPY.structure.state.LIDF.A = 1;
    spac.CANOPY.structure.state.LIDF.B = 1;
    EmeraldLand.CanopyOptics.inclination_angles!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;

    # By default, we use VerhoefLIDF method to compute LIDF, but we also support the use of Beta function.
    # To use the BetaLIDF, you need to change the parameter to BetaLIDF first.
    spac.CANOPY.structure.state.LIDF = EmeraldLand.Namespace.BetaLIDF{FT}();
    EmeraldLand.CanopyOptics.inclination_angles!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;
end;

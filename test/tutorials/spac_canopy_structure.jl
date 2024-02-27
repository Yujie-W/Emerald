using Emerald;
using Test;


@testset "Modify Canopy Structural Parameters" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_states!(config, spac);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.Namespace.t_aux!(config, spac);
    EmeraldLand.Namespace.s_aux!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    # Changing canopy structure may result in changes in other parameters, such as Vcmax profile.
    # Thus, it is not recommended to modify those parametes manually unless otherwise told to by developers.
    # It is highly recommended to use our embedded function prescribe_traits! to modify canopy structure.
    # Currently, function prescribe_traits! supports the modification of leaf area index and clumping index.
    EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 3, ci = 0.8);
    @test true;

    # Leaf inclination angle distribution is stored as a vector p_incl in field canopy.
    # By default, p_incl is a uniform distribution from 0° to 90° per 10°.
    # To change p_incl, we use the VerhoefLIDF model to update the distribution function.
    # For example, (A,B) =
    #     (0,0) gives uniform distribution,
    #     (-0.35,-0.15) gives spherical distribution,
    #     (-1,0) gives erectophile distribution,
    #     (1,0) gives planophile distribution,
    #     (0,-1) gives plagiophile distribution, and
    #     (0,1) gives extremophile distribution.
    spac.canopy.structure.trait.lidf.A = -0.35;
    spac.canopy.structure.trait.lidf.B = -0.15;
    EmeraldLand.Namespace.t_aux!(config, spac.canopy);
    EmeraldLand.Namespace.s_aux!(config, spac.canopy);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;

    # We also support beta function distribution. For example (A,B) =
    #     (1.000,1.000) gives uniform distribution,
    #     (1.930,1.101) gives spherical distribution,
    #     (2.770,1.172) gives erectophile distribution,
    #     (1.172,2.770) gives planophile distribution,
    #     (3.326,3.326) gives plagiophile distribution, and
    #     (0.433,0.433) gives extremophile distribution.
    spac.canopy.structure.trait.lidf = EmeraldLand.Namespace.BetaLIDF{FT}();
    spac.canopy.structure.trait.lidf.A = 1;
    spac.canopy.structure.trait.lidf.B = 1;
    EmeraldLand.Namespace.t_aux!(config, spac.canopy);
    EmeraldLand.Namespace.s_aux!(config, spac.canopy);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;

    # By default, we use VerhoefLIDF method to compute LIDF, but we also support the use of Beta function.
    # To use the BetaLIDF, you need to change the parameter to BetaLIDF first.
    spac.canopy.structure.trait.lidf = EmeraldLand.Namespace.BetaLIDF{FT}();
    EmeraldLand.Namespace.t_aux!(config, spac.canopy);
    EmeraldLand.Namespace.s_aux!(config, spac.canopy);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;
end;

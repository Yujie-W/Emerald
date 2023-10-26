using Emerald;
using Test;


@testset "Leaf Quantum Yields" begin
    FT = Float64;

    # Photosynthesis module can be used to compute leaf quantum yields
    # To do so, users can define and use the variables at photosystem level to avoid unnecessary calculations
    #     - Photosystem (C3VJP, C4VJP, or C3Cyto)
    #     - Air layer
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    ps = EmeraldLand.Namespace.LeafPhotosystem{FT}();
    air = EmeraldLand.Namespace.AirLayer{FT}();

    # Then run the following steps to compute quantum yields
    # These variables will be required as inputs
    #     - Leaf temperature
    #     - PAR that goes to photosystems (not the total PAR or APAR)
    #     - Internal CO₂ concentration [`Pa`]
    t_leaf = 298.15;
    ppar = 100.0;
    p_co2 = 20.0;
    EmeraldLand.Photosynthesis.photosystem_temperature_dependence!(ps, air, t_leaf);
    EmeraldLand.Photosynthesis.photosystem_electron_transport!(ps, ppar, p_co2);
    EmeraldLand.Photosynthesis.rubisco_limited_rate!(ps, p_co2);
    EmeraldLand.Photosynthesis.light_limited_rate!(ps);
    EmeraldLand.Photosynthesis.product_limited_rate!(ps, p_co2);
    EmeraldLand.Photosynthesis.colimit_photosynthesis!(ps);
    EmeraldLand.Photosynthesis.photosystem_coefficients!(ps, ppar);

    # Then the quantum yields can be read from the structure
    @test ps.auxil.ϕ_f > 0;
    @test ps.auxil.ϕ_p > 0;

    # If you want to change the physiological parameters, you can do so by changing the state variables (see the structure definition for more details), e.g.,
    #     - Vcmax25
    #     - Jmax25 (C3VJP only)
    #     - b₆f (C3Cyto only)
    #     - Vpmax25 (C4VJP only)
    #     - Rd25
    #     - Sustained NPQ (k_npq_sus; C3VJP and C4VJP only)
    # Note to change the memory temperature to a different value so that the temperature dependence
    # Afterward, you need to run the above steps again to update the quantum yields
    ps.state.v_cmax25 = 80;
    ps.state.j_max25 = 160;
    ps.state.r_d25 = 2;
    ps.auxil._t = 0;

    EmeraldLand.Photosynthesis.photosystem_temperature_dependence!(ps, air, t_leaf);
    EmeraldLand.Photosynthesis.photosystem_electron_transport!(ps, ppar, p_co2);
    EmeraldLand.Photosynthesis.rubisco_limited_rate!(ps, p_co2);
    EmeraldLand.Photosynthesis.light_limited_rate!(ps);
    EmeraldLand.Photosynthesis.product_limited_rate!(ps, p_co2);
    EmeraldLand.Photosynthesis.colimit_photosynthesis!(ps);
    EmeraldLand.Photosynthesis.photosystem_coefficients!(ps, ppar);

    @test ps.auxil.ϕ_f > 0;
    @test ps.auxil.ϕ_p > 0;
end;

using Emerald;
using Test;


@testset "Run SPAC Model" begin
    FT = Float64;

    # To run SPAC model, users need to define the following to proceed, and then initialize and run the model:
    #     - Configuration struct
    #     - SPAC struct
    # For details about how to modify the SPAC struct or configurations, please refer to other tutorials.
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_states!(config, spac);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.CanopyOptics.t_aux!(config, spac.canopy.structure.trait, spac.canopy.structure.t_aux);
    EmeraldLand.CanopyOptics.s_aux!(config, spac.canopy.structure.trait, spac.canopy.structure.t_aux, spac.canopy.sun_geometry.state, spac.canopy.sun_geometry.s_aux);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;
end;

using Emerald;
using Test;


@testset "Run SPAC Model" begin
    FT = Float64;

    # To run SPAC model, users need to define the following to proceed, and then initialize and run the model:
    #     - Configuration struct
    #     - SPAC struct
    # For details about how to modify the SPAC struct or configurations, please refer to other tutorials.
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;

    # change the LAI to 6 and SZA to 88.9°
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 6, sai = 0.6);
    spac.canopy.sun_geometry.state.sza = 88.9;
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;
end;

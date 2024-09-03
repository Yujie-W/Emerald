using Emerald;
using Test;


@testset "Kill a plant" begin
    FT = Float64;

    # To run SPAC model, users need to define the following to proceed, and then initialize and run the model:
    #     - Configuration struct
    #     - SPAC struct
    # For details about how to modify the SPAC struct or configurations, please refer to other tutorials.
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @test true;

    # Change carbon pool to negative and kill the plant
    spac.plant.pool.c_pool = -100;
    EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @test true;
end;

using Emerald;
using Test;


@testset "Run SPAC Model with binned PPAR" begin
    FT = Float64;

    # To run SPAC model, users need to define the following to proceed, and then initialize and run the model:
    #     - Configuration struct
    #     - SPAC struct
    # For details about how to modify the SPAC struct or configurations, please refer to other tutorials.
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, 3600);
    @test true;

    # bin PPAR to 10 and 20 bins
    config_10 = EmeraldLand.Namespace.SPACConfiguration(FT);
    config_10.DIM_PPAR_BINS = 10;
    config_20 = EmeraldLand.Namespace.SPACConfiguration(FT);
    config_20.DIM_PPAR_BINS = 20;
    spac_10 = EmeraldLand.Namespace.BulkSPAC(config_10);
    spac_20 = EmeraldLand.Namespace.BulkSPAC(config_20);
    EmeraldLand.SPAC.initialize_spac!(config_10, spac_10);
    EmeraldLand.SPAC.initialize_spac!(config_20, spac_20);
    EmeraldLand.SPAC.spac!(config_10, spac_10, 3600);
    EmeraldLand.SPAC.spac!(config_20, spac_20, 3600);
    @test true;

    # @info "Check the results" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.GPP(spac_10) EmeraldLand.SPAC.GPP(spac_20);
    # @info "Check the results" EmeraldLand.SPAC.TROPOMI_SIF740(config, spac) EmeraldLand.SPAC.TROPOMI_SIF740(config_10, spac_10) EmeraldLand.SPAC.TROPOMI_SIF740(config_20, spac_20);
end;

using Emerald;
using Test;


@testset "Change Stomatal Model" begin
    FT = Float64;

    # Emerald supports a number of stomatal models from empirical to optimality models. By default, the model is set to the Wang et al. (2020) and (2021) models.
    # Users can change the model to other models by setting the stomatal model in the SPAC struct.
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    for leaf in spac.plant.leaves
        leaf.flux.trait.stomatal_model = EmeraldLand.Namespace.MedlynSM{FT}();
    end;
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;
end;

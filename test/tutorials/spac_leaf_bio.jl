using Emerald;
using Test;


@testset "Modify Leaf Biophysical Parameters" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    # There are two ways to change leaf biophysical parameters.
    # The first method is to do it manually.
    # However, users need to run function leaf_spectra! to force the model update leaf level spectra.
    for leaf in spac.plant.leaves
        leaf.bio.trait.cab = 40;        # chlorophyll a and b content
        leaf.bio.trait.car = 10;        # carotenoid content
        leaf.bio.trait.cbc = 0.011;     # strutural carbon content
        leaf.bio.trait.lma = 0.012;     # leaf mass per area lma = cbc + pro
        leaf.bio.trait.pro = 0.001;     # protein content
    end;
    EmeraldLand.LeafOptics.plant_leaf_spectra!(config, spac);
    @test true;

    # The second method is to use the embedded function prescribe_traits! from submodule SPAC.
    # In this case, function leaf_spectra! runs automatically.
    # Note that it is recommended to change both at the same time, otherwise leaf_spectra! will run twice.
    # As to the supported options of the function prescribe_traits!, please check out the documentation.
    EmeraldLand.SPAC.prescribe_traits!(config, spac; cab = 40, car = 10);
    EmeraldLand.SPAC.dull_aux!(config, spac);
    @test true;
end;

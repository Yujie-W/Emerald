using Emerald;
using Test;


@testset "Modify Leaf Biophysical Parameters" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.MultiLayerSPAC(config);
    EmeraldLand.SPAC.initialize!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));

    # There are two ways to change leaf biophysical parameters.
    # The first method is to do it manually.
    # However, users need to run function leaf_spectra! to force the model update leaf level spectra.
    # Also, remember to change the cache variable _v_storage to a different value.
    # Otherwise, leaf spectra would not be updated.
    for leaf in spac.LEAVES
        leaf.NS.bio.state.cab = 40;        # chlorophyll a and b content
        leaf.NS.bio.state.car = 10;        # carotenoid content
        leaf.NS.bio.state.cbc = 0.011;     # strutural carbon content
        leaf.NS.bio.state.lma = 0.012;     # leaf mass per area lma = cbc + pro
        leaf.NS.bio.state.pro = 0.001;     # protein content
    end;
    EmeraldLand.LeafOptics.leaf_spectra!(config, spac);
    @test true;

    # The second method is to use the embedded function update! from submodule SPAC.
    # In this case, function leaf_spectra! runs automatically.
    # Note that it is recommended to change both at the same time, otherwise leaf_spectra! will run twice.
    # As to the supported options of the function update!, please check out the documentation.
    EmeraldLand.SPAC.update!(config, spac; cab = 40, car = 10);
    @test true;
end

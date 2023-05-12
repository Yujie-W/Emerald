@testset "Modify Leaf Biophysical Parameters" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.MultiLayerSPAC{FT}();
    EmeraldLand.SPAC.initialize!(spac, config);
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    # There are two ways to change leaf biophysical parameters.
    # The first method is to do it manually.
    # However, users need to run function leaf_spectra! to force the model update leaf level spectra.
    # Also, remember to change the cache variable _v_storage to a different value.
    # Otherwise, leaf spectra would not be updated.
    for leaf in spac.LEAVES
        leaf.BIO._v_storage = 0;    # change _v_storage to 0 so as to trigger function leaf_spectra!
        leaf.BIO.cab = 40;          # chlorophyll a and b content
        leaf.BIO.car = 10;          # carotenoid content
        leaf.BIO.cbc = 0.011;       # strutural carbon content
        leaf.BIO.lma = 0.012;       # leaf mass per area lma = cbc + pro
        leaf.BIO.pro = 0.001;       # protein content
    end;
    EmeraldLand.LeafOptics.leaf_spectra!(spac, config);
    @test true;

    # The second method is to use the embedded function update! from submodule SPAC.
    # In this case, function leaf_spectra! runs automatically.
    # Note that it is recommended to change both at the same time, otherwise leaf_spectra! will run twice.
    # As to the supported options of the function update!, please check out the documentation.
    EmeraldLand.SPAC.update!(spac, config; cab = 40, car = 10);
    @test true;
end

@testset "Run Leaf Level Spectra" begin
    FT = Float64;

    # LeafOptics package can be used in stand-alone mode without loading unnecessary modules.
    # To do so, users need to define and use the variables at leaf level:
    #     - Leaf biophysical parameters such as chlorophyll
    #     - Wavelength parameter set that define the wavelength bins
    #     - Hyperspectral absorption feature of leaf biophysical traits
    #     - Leaf water content in [mol m⁻²] (yes, as an input)
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    bio = EmeraldLand.Namespace.HyperLeafBio(config);
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;

    # Note that, to avoid unnecessary computing, leaf water content is saved in bio._v_storage.
    # If you change leaf biophysical traits such as chlorophyll content, you need to set _v_storage to 0 (or any different water content value).
    bio.state.cab = 20;
    bio.state.car = 5;
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;
end

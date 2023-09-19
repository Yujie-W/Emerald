@testset "Run Leaf Level Spectra" begin
    FT = Float64;

    # LeafOptics package can be used in stand-alone mode without loading unnecessary modules.
    # To do so, users need to define and use the variables at leaf level:
    #     - Leaf biophysical parameters such as chlorophyll
    #     - Wavelength parameter set that define the wavelength bins
    #     - Hyperspectral absorption feature of leaf biophysical traits
    #     - Leaf water content in [mol m⁻²] (yes, as an input)
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    bio = EmeraldLand.Namespace.HyperspectralLeafBiophysics(config);
    wls = EmeraldLand.Namespace.WaveLengthSet{FT}();
    spectra = EmeraldLand.Namespace.ReferenceSpectra{FT}();
    EmeraldLand.LeafOptics.leaf_spectra!(bio, wls, spectra, FT(5));
    @test true;

    # Note that, to avoid unnecessary computing, leaf water content is saved in bio._v_storage.
    # If you change leaf biophysical traits such as chlorophyll content, you need to set _v_storage to 0 (or any different water content value).
    bio.cab = 20;
    bio.car = 5;
    bio._v_storage = 0;
    EmeraldLand.LeafOptics.leaf_spectra!(bio, wls, spectra, FT(5));
    @test true;

    # Users can also emulate broadband simulations in hyperspectral mode by providing
    #     - PAR reflectance
    #     - NIR reflectance
    #     - PAR transmittance
    #     - NIR transmittance
    EmeraldLand.LeafOptics.leaf_spectra!(bio, wls, FT(0.1), FT(0.2), FT(0.1), FT(0.2));
    @test true;
end

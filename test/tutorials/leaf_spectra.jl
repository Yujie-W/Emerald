using Emerald;
using Test;


@testset "Leaf Level Spectra" begin
    FT = Float64;

    # LeafOptics package can be used in stand-alone mode without loading unnecessary modules.
    # To do so, users need to define and use the variables at leaf level:
    #     - Leaf biophysical parameters such as chlorophyll
    #     - Wavelength parameter set that define the wavelength bins
    #     - Hyperspectral absorption feature of leaf biophysical traits
    #     - Leaf water content in [mol m⁻²] (yes, as an input)
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    bio = EmeraldLand.Namespace.LeafBio(config);
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;

    # Change the leaf biophysical parameters
    bio.trait.cab = 20;
    bio.trait.car = 5;
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;
end;

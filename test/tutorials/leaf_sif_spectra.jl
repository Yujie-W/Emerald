using Emerald;
using Test;


@testset "Examine SIF Spectra" begin
    FT = Float64;

    # Use LeafOptics module to calculate SIF spectra at different conditions
    # Here, we use 1 nm wavelength resolution for better accuracy
    # To do so, users need to define and use the variables at leaf level:
    #     - General SPAC configuration
    #     - Leaf biophysical parameters such as chlorophyll
    #     - Leaf water content in [mol m⁻²] (yes, as an input)
    # Note that in the config struct, there are three features that affect the SIF calculation
    #     Φ_SIF_CUTOFF      Int     how to cut off the SIF emission spectrum at the excitation wavelength (default 0)
    #                               0: no cut off
    #                               1: sharp cut off
    #                               2: sigmoid used in SCOPE
    #     Φ_SIF_RESCALE     Bool    whether to rescale the SIF emission PDF after cut off (default true)
    #     Φ_SIF_WL          Bool    whether to partition the SIF emission PDF based on the excitation wavelength (default true)
    config = EmeraldLand.Namespace.SPACConfiguration{FT}(DATASET = EmeraldLand.Namespace.LAND_2021_1NM, Φ_SIF_CUTOFF = 0, Φ_SIF_RESCALE = true, Φ_SIF_WL = true);
    bio = EmeraldLand.Namespace.LeafBio(config);
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;

    # Now the SIF spectra can be computed using the stored matrices mat_b and mat_f for any given incoming radiation spectrum
    # Here we define a blue light and red light (only in the SIF excitation spectrum)
    rad_b = zeros(FT, length(config.SPECTRA.IΛ_SIFE));
    rad_r = zeros(FT, length(config.SPECTRA.IΛ_SIFE));
    rad_b[450 .<= config.SPECTRA.Λ_SIFE .<= 480] .= 1;
    rad_r[620 .<= config.SPECTRA.Λ_SIFE .<= 650] .= 1;
    sif_b_b = bio.auxil.mat_b * rad_b;
    sif_f_b = bio.auxil.mat_f * rad_b;
    sif_b_r = bio.auxil.mat_b * rad_r;
    sif_f_r = bio.auxil.mat_f * rad_r;
    @test true;

    # If you change leaf biophysical traits such as chlorophyll content, you can change them in the state field of the bio
    bio.trait.cab = 20;
    bio.trait.car = 5;
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;
end;

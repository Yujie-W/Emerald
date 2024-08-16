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
    config = EmeraldLand.Namespace.SPACConfiguration(FT; dataset = EmeraldLand.Namespace.OLD_PHI_2021_1NM);
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
    #=
    using Emerald;
    using DataFrames;
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration(FT; dataset = EmeraldLand.Namespace.OLD_PHI_2021_1NM);
    config.Φ_SIF_CUTOFF = 2;
    bio = EmeraldLand.Namespace.LeafBio(config);
    bio.trait.cab = 15;
    bio.trait.car = 3;
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    rad_b = zeros(FT, length(config.SPECTRA.IΛ_SIFE));
    rad_r = zeros(FT, length(config.SPECTRA.IΛ_SIFE));
    rad_b[450 .<= config.SPECTRA.Λ_SIFE .<= 480] .= 1;
    rad_r[620 .<= config.SPECTRA.Λ_SIFE .<= 650] .= 1;
    sif_b_b = bio.auxil.mat_b * rad_b;
    sif_f_b = bio.auxil.mat_f * rad_b;
    sif_b_r = bio.auxil.mat_b * rad_r;
    sif_f_r = bio.auxil.mat_f * rad_r;
    df = DataFrame(WL = config.SPECTRA.Λ_SIF, SIF_BW_B = sif_b_b, SIF_FW_B = sif_f_b, SIF_BW_R = sif_b_r, SIF_FW_R = sif_f_r);
    EmeraldIO.Text.save_csv!("sif-shapes.csv", df);
    =#
    @test true;

    # If you change leaf biophysical traits such as chlorophyll content, you can change them in the state field of the bio
    bio.trait.cab = 20;
    bio.trait.car = 5;
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;

    # the default SIF model is my new SIF model, but you can change it to the old model in Fluspect-B
    config = EmeraldLand.Namespace.SPACConfiguration(FT; dataset = EmeraldLand.Namespace.OLD_PHI_2021_1NM);
    bio = EmeraldLand.Namespace.LeafBio(config);
    bio.trait.SIF_METHOD = EmeraldLand.Namespace.SIFMatrixDoublingMethod();
    EmeraldLand.LeafOptics.leaf_spectra!(config, bio, FT(5));
    @test true;
end;

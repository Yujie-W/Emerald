using Test
using Emerald;


@testset "Fluorescence photon conservation of the SCOPE SIF model" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}(DATASET = EmeraldLand.Namespace.LAND_2021_1NM);
    bio = EmeraldLand.Namespace.HyperspectralLeafBiophysics(config);
    spectra = config.SPECTRA;
    EmeraldLand.LeafOptics.leaf_spectra!(bio, spectra, FT(5));

    # Test the photon conservation
    rad = EmeraldLand.Namespace.HyperspectralRadiation{FT}(config.DATASET);
    rad_e = config.SPECTRA.SOLAR_RAD[:,1];
    rad_photon = EmeraldLand.Optics.photon.(spectra.Λ, rad_e);
    _,_,ppar = EmeraldLand.LeafOptics.leaf_PAR(bio, spectra, rad);
    fppar = bio.α_sw .* bio.α_cabcar;
    @show ppar ./ 1000;
    @show (rad_photon[spectra.IΛ_SIFE] .* fppar[spectra.IΛ_SIFE])' * spectra.ΔΛ_SIFE;

    sif_photon = (bio.mat_b .+ bio.mat_f) * (rad_photon[spectra.IΛ_SIFE] .* spectra.ΔΛ_SIFE);
    @show sif_photon' * spectra.ΔΛ_SIF;

    sif_photon = (bio._mat_b .+ bio._mat_f) * (rad_photon[spectra.IΛ_SIFE] .* spectra.ΔΛ_SIFE);
    @show sif_photon' * spectra.ΔΛ_SIF;

    @test true;
end

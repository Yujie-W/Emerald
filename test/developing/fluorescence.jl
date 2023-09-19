@testset "Fluorescence photon conservation" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}(DATASET = EmeraldLand.Namespace.LAND_2021_1NM);
    bio = EmeraldLand.Namespace.HyperspectralLeafBiophysics(config);
    wls = config.WLSET;
    lha = config.LHA;
    EmeraldLand.LeafOptics.leaf_spectra!(bio, wls, lha, FT(5));





    # Test the photon conservation
    rad = config.RAD_SW_REF;
    rad_e = rad.e_direct .+ rad.e_diffuse;
    rad_photon = EmeraldLand.Optics.photon.(wls.Λ, rad_e);
    _,_,ppar = EmeraldLand.LeafOptics.leaf_PAR(bio, wls, rad);
    fppar = bio.α_sw .* bio.α_cabcar;
    @show ppar ./ 1000;
    @show (rad_photon[wls.IΛ_SIFE] .* fppar[wls.IΛ_SIFE])' * wls.ΔΛ_SIFE;

    sif_photon = (bio.mat_b .+ bio.mat_f) * (rad_photon[wls.IΛ_SIFE] .* wls.ΔΛ_SIFE);
    @show sif_photon' * wls.ΔΛ_SIF;

    sif_photon = (bio._mat_b .+ bio._mat_f) * (rad_photon[wls.IΛ_SIFE] .* wls.ΔΛ_SIFE);
    @show sif_photon' * wls.ΔΛ_SIF;







    @test true;
end

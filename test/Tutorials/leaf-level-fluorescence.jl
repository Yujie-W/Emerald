using PkgUtility.UniversalConstants: energy_to_photon
using Test

import Emerald.LeafOptics as ELO
import Emerald.Namespace as ENS
import Emerald.Photosynthesis as EPS
import Emerald.ResearchTools as ERT


@testset "Leaf Level Fluorescence" verbose = true begin
    @testset "C3 - Jmax + Platespect" verbose = true begin
        # these config settings are the default methods, I am just being explicit here
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C3_AC_METHOD = ENS.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = ENS.AjMethodC3JmaxPi();
        config.METHODS.C3_AP_METHOD = ENS.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = ENS.ColimitJCLM(Float64);
        config.METHODS.FLUORESCENCE_METHOD_C3 = ENS.KNFluorescenceModel{Float64}();
        cache = ERT.LeafLevelSetup.leaf_level_spac_cache(config);
        leaf = ERT.LeafLevelSetup.leaf_level_leaf(config, "C3");
        air = ENS.AirLayer{Float64}();

        # initialize the leaf water storage to full manually, otherwise, you may customize it as needed
        leaf.capacitor.state.v_storage = leaf.capacitor.trait.v_max;
        ELO.leaf_spectra!(config, leaf.bio, cache, leaf.capacitor.state.v_storage);
        @test all(0 .< leaf.bio.auxil.ρ_leaf .< 1);
        @test all(0 .< leaf.bio.auxil.τ_leaf .< 1);
        @test all(0 .< leaf.bio.auxil.mat_b .< 1);
        @test all(0 .< leaf.bio.auxil.mat_f .< 1);

        # compute backward and forward fluorescence using default radiation
        rad = config.CONSTANTS.SPECTRA.SOLAR_RAD[:,1];
        rad_excite = rad[config.CONSTANTS.SPECTRA.IΛ_SIFE,1];
        rad_ppar = rad .* leaf.bio.auxil.α_leaf .* leaf.bio.auxil.f_ppar;
        photon_ppar = energy_to_photon.(config.CONSTANTS.SPECTRA.Λ, rad_ppar .* 1e-3);

        # compute the photosynthesis and fluorescence yield (leaf.bio.auxil.f_ppar is the same as leaf.bio.auxil.f_sife by default)
        p_i = 20.0;                                                 # Pa
        ppar = photon_ppar' * config.CONSTANTS.SPECTRA.ΔΛ * 1e6;    # μmol m⁻² s⁻¹
        t = 298.15;                                                 # K
        EPS.photosynthesis!(config, cache, leaf.photosystem, air, [p_i,], [ppar,], t);

        # compare the fluorescence outputs
        sif_b = (leaf.bio.auxil.mat_b * rad_excite) .* leaf.photosystem.auxil.ϕ_f;
        sif_f = (leaf.bio.auxil.mat_f * rad_excite) .* leaf.photosystem.auxil.ϕ_f;
        photon_b = energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIF, sif_b .* 1e-3);
        photon_f = energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIF, sif_f .* 1e-3);
        sif_photon_leaf = (photon_b + photon_f)' * config.CONSTANTS.SPECTRA.ΔΛ_SIF * 1e6;  # μmol m⁻² s⁻¹
        sif_photon_chl = ppar * leaf.photosystem.auxil.ϕ_f[1];                             # μmol m⁻² s⁻¹
        @assert sif_photon_chl > sif_photon_leaf > 0;
    end;
end;

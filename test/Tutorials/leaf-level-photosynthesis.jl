using Test

import Emerald.Namespace as ENS
import Emerald.Photosynthesis as EPS
import Emerald.ResearchTools as ERT


@testset "Leaf Level Photosynthesis" verbose = true begin
    @testset "C3 - Jmax" begin
        # these config settings are the default methods, I am just being explicit here
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C3_AC_METHOD = ENS.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = ENS.AjMethodC3JmaxPi();
        config.METHODS.C3_AP_METHOD = ENS.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = ENS.ColimitJCLM(Float64);
        config.METHODS.FLUORESCENCE_METHOD = ENS.KNFluorescenceModel{Float64}();
        cache = ERT.LeafLevelSetup.leaf_level_spac_cache(config);
        lps = ERT.LeafLevelSetup.leaf_level_photosystem(Float64, "C3");
        air = ENS.AirLayer{Float64}();
        p_i = 20.0;     # Pa
        ppar = 1000.0;  # μmol m⁻² s⁻¹
        t = 298.15;     # K
        EPS.photosynthesis!(config, cache, lps, air, [p_i,], [ppar,], t);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C3 - Vqmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.C3_AC_METHOD = ENS.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = ENS.AjMethodC3VqmaxPi();
        config.METHODS.C3_AP_METHOD = ENS.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = ENS.SerialColimit();
        config.METHODS.FLUORESCENCE_METHOD = ENS.CytochromeFluorescenceModel();
        cache = ERT.LeafLevelSetup.leaf_level_spac_cache(config);
        lps = ERT.LeafLevelSetup.leaf_level_photosystem(Float64, "C3");
        air = ENS.AirLayer{Float64}();
        p_i = 20.0;     # Pa
        ppar = 1000.0;  # μmol m⁻² s⁻¹
        t = 298.15;     # K
        EPS.photosynthesis!(config, cache, lps, air, [p_i,], [ppar,], t);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;
end;

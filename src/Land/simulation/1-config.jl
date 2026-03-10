"""

    site_config(settings::Union{Dict,OrderedDict})

Create a SPAC configuration struct, given
- `gmd` Dictionary of GriddingMachine data in a grid

"""
function site_config(settings::Union{Dict,OrderedDict})
    config = SPACConfig(settings["FT"]);
    config.CONFIG_INFO.MESSAGE_LEVEL = settings["MESSAGE_LEVEL"];

    # set up the default features
    config.FEATURES.ALLOW_LEAF_REGROWTH = false;
    config.FEATURES.ALLOW_LEAF_SHEDDING = false;
    config.FEATURES.ALLOW_XYLEM_GROWTH = false;
    config.FEATURES.EFFECTIVE_LEAF_SPECTRA = false;
    config.FEATURES.ENABLE_DROUGHT_LEGACY = false;
    config.FEATURES.ENABLE_REF = true;
    config.FEATURES.ENABLE_SIF = true;
    config.FEATURES.UNLIMITED_NSC_POOL = true;

    # set up the default methods
    config.METHODS.STOMATAL_MODEL = Namespace.WangSM{settings["FT"]}();

    # set up the photosynthesis model
    if settings["C3_MODEL"] == "Jmax"
        config.METHODS.C3_AC_METHOD = Namespace.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = Namespace.AjMethodC3JmaxPi();
        config.METHODS.C3_AP_METHOD = Namespace.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = Namespace.ColimitJCLM(settings["FT"]);
        # TODO: these results are from the unpublished research I am working on, will make them official parameter sets when the paper is accepted
        config.METHODS.TD_VCMAX_C3.ΔHA = 63000;
        config.METHODS.TD_VCMAX_C3.ΔHD = 204000;
        config.METHODS.TD_JMAX.ΔHA = 50000;
        config.METHODS.TD_JMAX.ΔHD = 201000;
        config.METHODS.TD_Γ.VAL_REF = 4.67;
        config.METHODS.TD_Γ.ΔHA = 11800;
    elseif settings["C3_MODEL"] == "Vqmax"
        config.METHODS.C3_AC_METHOD = Namespace.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = Namespace.AjMethodC3VqmaxPi();
        config.METHODS.C3_AP_METHOD = Namespace.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = Namespace.SerialColimit();
        # TODO: these results are from the unpublished research I am working on, will make them official parameter sets when the paper is accepted
        config.METHODS.TD_VCMAX_C3.ΔHA = 63000;
        config.METHODS.TD_VCMAX_C3.ΔHD = 204000;
        config.METHODS.TD_KQ = Namespace.ArrheniusPeak{Float64}(298.15, 300, 28500, 223500, 700);
        config.METHODS.TD_Γ.VAL_REF = 4.56;
        config.METHODS.TD_Γ.ΔHA = 11800;
        config.METHODS.TD_ηC.VAL_REF = 4 * 3 / 14;
        config.METHODS.TD_ηC.ΔHA = 28500;
        config.METHODS.TD_ηC.ΔHD = 223500;
        config.METHODS.TD_ηC.ΔSV = 700;
        config.METHODS.TD_ηL.VAL_REF = 3 * 3 / 14;
        config.METHODS.TD_ηL.ΔHA = 28500;
        config.METHODS.TD_ηL.ΔHD = 223500;
        config.METHODS.TD_ηL.ΔSV = 700;
    else
        pretty_display!("C3 photosynthesis model not recognized: $(settings["C3_MODEL"]), use testing setting instead...", "twarn");
    end;

    # set up the C4 photosynthesis model
    if settings["C4_MODEL"] == "Vcmax"
        config.METHODS.C4_AC_METHOD = Namespace.AcMethodC4Vcmax();
        config.METHODS.C4_AJ_METHOD = Namespace.AjMethodC4JPSII();
        config.METHODS.C4_AP_METHOD = Namespace.ApMethodC4VcmaxPi();
        config.METHODS.TD_R_C4 = Namespace.RespirationTDCLMC4(settings["FT"]);
        config.METHODS.TD_VCMAX_C4 = Namespace.VcmaxTDCLMC4(settings["FT"]);
    elseif settings["C4_MODEL"] == "Vpmax"
        config.METHODS.C4_AC_METHOD = Namespace.AcMethodC4Vcmax();
        config.METHODS.C4_AJ_METHOD = Namespace.AjMethodC4JPSII();
        config.METHODS.C4_AP_METHOD = Namespace.ApMethodC4VpmaxPi();
        config.METHODS.TD_R_C4 = Namespace.RespirationTDCLMC4(settings["FT"]);
        config.METHODS.TD_VCMAX_C4 = Namespace.VcmaxTDCLMC4(settings["FT"]);
        config.METHODS.TD_VPMAX = Namespace.VpmaxTDBoyd(settings["FT"]);
    end;

    # set up the fluorescence model for C3 plants
    config.METHODS.FLUORESCENCE_METHOD_C3 = if settings["C3_ΦF_MODEL"] == "KN"
        Namespace.KNFluorescenceModel{settings["FT"]}()
    elseif settings["C3_ΦF_MODEL"] == "QL"
        Namespace.QLFluorescenceModelHanC3(settings["FT"])
    elseif settings["C3_ΦF_MODEL"] == "B6F"
        @assert settings["C3_MODEL"] == "Vqmax" "The J3B fluorescence model can only be used along with the Vqmax photosynthesis model";
        Namespace.CytochromeFluorescenceModel()
    else
        error("When C3_MODEL is FvCB, C3_ΦF_MODEL must be either KN, QL, or B6F, but got $(settings["C3_ΦF_MODEL"])...")
    end;

    # set up the fluorescence model for C4 plants
    config.METHODS.FLUORESCENCE_METHOD_C4 = if settings["C4_ΦF_MODEL"] == "KN"
        Namespace.KNFluorescenceModel{settings["FT"]}()
    elseif settings["C4_ΦF_MODEL"] == "QL"
        Namespace.QLFluorescenceModelHanC4(settings["FT"])
    else
        error("When C4_MODEL is Vcmax, C4_ΦF_MODEL must be either KN or QL, but got $(settings["C4_ΦF_MODEL"])...")
    end;

    # set up soil albedo model
    if settings["SOIL_ALBEDO_MODEL"] == "HyperspectralCliMA"
        config.METHODS.SOIL_ALBEDO = Namespace.SoilAlbedoHyperspectralCLIMA();
    elseif settings["SOIL_ALBEDO_MODEL"] == "HyperspectralAsh"
        config.METHODS.SOIL_ALBEDO = Namespace.SoilAlbedoHyperspectralAsh();
        # TODO add the data for soil albedo from the ASD experiments
    end;

    return config
end;

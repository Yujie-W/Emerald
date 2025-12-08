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
    if settings["C3_MODEL"] == "FvCB"
        config.METHODS.C3_AC_METHOD = Namespace.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = Namespace.AjMethodC3JmaxPi();
        config.METHODS.C3_AP_METHOD = Namespace.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = Namespace.ColimitJCLM(settings["FT"]);
    elseif settings["C3_MODEL"] == "J3B"
        config.METHODS.C3_AC_METHOD = Namespace.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = Namespace.AjMethodC3VqmaxPi();
        config.METHODS.C3_AP_METHOD = Namespace.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = Namespace.SerialColimit{settings["FT"]}();
    else
        pretty_display!("C3 photosynthesis model not recognized: $(settings["C3_MODEL"]), use testing setting instead...", "twarn");
    end;

    return config
end;

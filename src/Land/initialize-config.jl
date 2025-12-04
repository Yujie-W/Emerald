"""

    site_config(gmd::Dict)

Create a SPAC configuration struct, given
- `gmd` Dictionary of GriddingMachine data in a grid

"""
function site_config(gmd::Dict)
    config = SPACConfig(gmd["FT"]);
    config.CONFIG_INFO.MESSAGE_LEVEL = gmd["MESSAGE_LEVEL"];

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
    config.METHODS.STOMATAL_MODEL = WangSM{gmd["FT"]}();

    return config
end;

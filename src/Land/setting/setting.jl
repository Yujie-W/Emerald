"""

    land_model_settings(; mode::String = "testing")

Create a dictionary of Emerald land settings, given
- `mode` Mode of the settings, e.g., "testing", "default", or "cytochrome"

"""
function land_model_settings(; mode::String = "testing")
    settings = OrderedDict{String,Any}(
        # Emerald version
        "EMERALD_VERSION"      => "b01",
        "CONFIG_TAG"           => "testing",

        # general settings
        "FT"                   => Float64,
        "NX"                   => 1,
        "GM_VERSION"           => "gm2",
        "WD_VERSION"           => "wd1",
        "MESSAGE_LEVEL"        => 0,
        "TIME_STEP"            => 3600,

        # SPAC settings
        "C3_MODEL"             => "FvCB",

        # threading settings (default is 75% of CPU cores)
        "GRID_THREADS"         => 40,
        "SIMU_THREADS"         => Int(ceil(Sys.CPU_THREADS * 0.75)),
        "REMOVE_WHEN_DONE"     => true,

        # saving settings related to the global NetCDF output files
        "VARIABLES_TO_SAVE"    => String["GPP", "ET", "PCI", "PPAR", "SIF740", "ΦF", "ΦP", "ΣSIF", "ΣSIF_CHL", "ΣSIF_LEAF"],
        "VARIABLES_TO_COMBINE" => String["GPP", "ET", "PCI", "PPAR", "SIF740", "ΦF", "ΦP", "ΣSIF", "ΣSIF_CHL", "ΣSIF_LEAF"],

        # testing settings (by default, run the model for 10 days in the middle of a year)
        "SIMULATION_PERIOD"    => 4321:4344,
    );

    # if mode is default
    if mode == "default"
        settings["CONFIG_TAG"] = "default";
        settings["SIMULATION_PERIOD"] = :;
    end;

    # if mode is cytochrome
    if mode == "cytochrome"
        settings["CONFIG_TAG"] = "cytochrome";
        settings["C3_MODEL"] = "J3B";
        settings["SIMULATION_PERIOD"] = :;
    end;

    return settings
end;

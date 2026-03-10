"""

    land_model_settings(; mode::String = "testing")

Create a dictionary of Emerald land settings, given
- `mode` Mode of the settings, e.g., "testing", "default", "cytochrome", "ash"

"""
function land_model_settings(; mode::String = "testing")
    # default settings
    settings = OrderedDict{String,Any}(
        # Emerald version
        "EMERALD_VERSION"   => "b01",
        "CONFIG_TAG"        => mode,

        # general settings
        "FT"                => Float64,
        "NX"                => 1,
        "GM_VERSION"        => "gm2",
        "WD_VERSION"        => "wd1",
        "MESSAGE_LEVEL"     => 0,
        "TIME_STEP"         => 3600,

        # SPAC settings
        "C3_MODEL"          => "FvCB",
        "C3_ΦF_MODEL"       => "KN",
        "MAX_LAI_LAYERING"  => false,
        "SOIL_ALBEDO_MODEL" => "HyperspectralCliMA",

        # threading settings (default is 75% of CPU cores)
        "GRID_THREADS"      => min(Int(ceil(Sys.CPU_THREADS * 0.75)), 40),
        "SIMU_THREADS"      => min(Int(ceil(Sys.CPU_THREADS * 0.75)), 480),
        "REMOVE_WHEN_DONE"  => true,

        # saving settings related to the global NetCDF output files
        "VARIABLES_TO_SAVE" => String["GPP", "ET", "SIF740"],

        # testing settings (by default, run the model for 10 days in the middle of a year)
        "SIMULATION_PERIOD" => :,
    );

    # if mode is default
    if mode == "default"
        return settings
    end;

    # if mode contains "testing"
    if occursin("testing", mode)
        settings["SIMULATION_PERIOD"] = 4321:4344;
    end;

    # if mode contains SIF (this is meant for Christian's SIF experiments)
    if occursin("SIF", mode)
        settings["MAX_LAI_LAYERING"] = true;
        for vn in ["PCI", "PPAR", "SIF740", "ΦF", "ΦP", "ΣSIF", "ΣSIF_CHL", "ΣSIF_LEAF"]
            if !(vn in settings["VARIABLES_TO_SAVE"])
                push!(settings["VARIABLES_TO_SAVE"], vn);
            end;
        end;

        return settings
    end;

    # if mode is cytochrome
    if occursin("cytochrome", mode)
        settings["C3_MODEL"] = "J3B";

        return settings
    end;

    # if mode is ash
    if occursin("ash", mode)
        settings["SOIL_ALBEDO_MODEL"] = "HyperspectralAsh";

        return settings
    end;

    # default is testing, and print warning if unrecognized mode
    if mode != "testing"
        @warn "Unrecognized mode '$mode'; use 'testing' settings...";
    end;

    return settings
end;

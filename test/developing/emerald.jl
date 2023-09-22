#
# Version B4 (default setting; run it on Curry)
#
using Emerald;
using ProgressMeter;

FT = Float64;

EMERALD_VER = "b4";
GMDATA_VER = "gm2";
WDRIVER_VER = "wd1";

if gethostname()[1:5] == "curry"
    # use 64 threads on curry
    EmeraldEarth.add_threads!(64, FT);

    # define the datasets, grids, states, and weather drivers
    dts = EmeraldEarth.LandDatasets{FT}(GMDATA_VER, 2020);
    mat = EmeraldEarth.gm_grids(dts);
    sts = Matrix{Union{Nothing,EmeraldLand.Namespace.MultiLayerSPACState{FT}}}(nothing, size(dts.t_lm));
    wdr = EmeraldData.ERA5.ERA5SingleLevelsDriver();

    # preload the weather driver
    wds = EmeraldEarth.weather_drivers(dts, wdr);

    # run global scale simulation for one day
    rm("global-test.nc"; force = true);
    @showprogress for ind in 1:24
        wdx = EmeraldEarth.wd_grids(dts, wds, ind; displaying = false);
        sts = EmeraldEarth.simulation!(mat, wdx, sts; displaying = false);
        EmeraldEarth.save_simulations!("global-test.nc", sts, ind; displaying = false);
    end;
else
    @warn "You are recommened to use Curry to run the global simulations!";
end;

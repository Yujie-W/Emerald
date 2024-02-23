#
# Version B5 (default setting; run it on Curry)
#
using Emerald;
using ProgressMeter;

FT = Float64;

EMERALD_VER = "b5";
GMDATA_VER = "gm3";
WDRIVER_VER = "wd1";

if gethostname()[1:5] == "curry"
    # define the datasets, grids, states, and weather drivers
    dts = EmeraldData.GlobalDatasets.LandDatasets{FT}(GMDATA_VER, 2019);
    mat = EmeraldData.GlobalDatasets.grid_dict_mat(dts);
    sts = Matrix{Union{Nothing,EmeraldLand.Namespace.BulkSPACStates{FT}}}(nothing, size(dts.t_lm));
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1);


    # use multiple threads on curry
    EmeraldEarth.add_threads!(32, FT);


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

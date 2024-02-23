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
    gms = EmeraldData.GlobalDatasets.grid_dict_mat(dts);
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1);
    iss = EmeraldData.WeatherDrivers.initial_soil_skin_states(WDRIVER_VER, 2019, 1, 1);

    # use multiple threads on curry
    EmeraldEarth.add_threads!(32, FT);
    wdi = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    mss = EmeraldData.WeatherDrivers.initial_states_snapshot(iss);

    EmeraldEarth.initial_states(gms, wdi, mss)


    EmeraldEarth.initial_state(gms[1], wdi[1], mss[1])



    sts = Matrix{Union{Nothing,EmeraldLand.Namespace.BulkSPACStates{FT}}}(nothing, size(dts.t_lm));
    wdi = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);




    # preload the weather driver
    wds = EmeraldEarth.weather_drivers(dts, wdr);

    # run global scale simulation for one day
    rm("global-test.nc"; force = true);
    @showprogress for ind in 1:24
        wdx = EmeraldEarth.wd_grids(dts, wds, ind; displaying = false);
        sts = EmeraldEarth.simulation!(gms, wdx, sts; displaying = false);
        EmeraldEarth.save_simulations!("global-test.nc", sts, ind; displaying = false);
    end;
else
    @warn "You are recommened to use Curry to run the global simulations!";
end;

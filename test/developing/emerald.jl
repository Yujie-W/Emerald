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
    # define the
    #     dts: land datasets
    #     gms: griddingmachine dataset matrix
    #     wds: preloaded weather drivers
    #     iss: initial soil and skin states
    #     wdi: weather drivers snapshot matrix
    #     mss: initial states snapshot matrix
    dts = EmeraldData.GlobalDatasets.LandDatasets{FT}(GMDATA_VER, 2019);
    gms = EmeraldData.GlobalDatasets.grid_dict_mat(dts; vegetation_only = true);
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, 1);
    iss = EmeraldData.WeatherDrivers.initial_soil_skin_states(WDRIVER_VER, 2019, 1, 1);
    wd1 = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    mss = EmeraldData.WeatherDrivers.initial_states_snapshot(iss);

    # use multiple threads on curry and initialize the states matrix with multiple threads and run the model for one time step
    EmeraldEarth.add_threads!(160, FT);
    ini_sts = EmeraldEarth.initial_states(gms, wd1, mss);
    sts_2nd = EmeraldEarth.global_simulations!(gms, wd1, ini_sts);

    # run the model for the second time step
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, 2);
    wd2 = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    sts_3rd = EmeraldEarth.global_simulations!(gms, wd2, sts_2nd);
    EmeraldEarth.global_simulations!(gms, wd2, sts_2nd; single_thread = true, single_thread_regions = (:,:));





    # run the global simulations for 10 time steps
    for i in 1:10
        @time wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, i);
        @time wdi = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
        @time new_sts = EmeraldEarth.global_simulations!(gms, wdi, new_sts);
    end;




    EmeraldEarth.global_simulations!(gms, wd1, ini_sts);
    EmeraldEarth.global_simulations!(gms, wd1, ini_sts; single_thread = true, single_thread_regions = (90:119,:));

    begin
        i = 97;
        j = 130;
        gmij = gms[i,j];
        wdij = wdi[i,j];
        ssij = mss[i,j];
        stij = deepcopy(sts[i,j]);
        EmeraldEarth.grid_simulation!(gmij, wdij, stij);
    end;





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

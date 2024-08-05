#
# Version B5 (default setting; run it on Curry)
#
using Emerald;
using ProgressMeter;

FT = Float64;

EMERALD_VER = "b5";
GMDATA_VER = "gm2";
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
    sts_1st = EmeraldEarth.initial_states(gms, wd1, mss);
    sts_2nd = EmeraldEarth.global_simulations!(gms, wd1, sts_1st);

    # run the model for the second time step
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, 2);
    wd2 = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    sts_3rd = EmeraldEarth.global_simulations!(gms, wd2, sts_2nd);

    # run the model for the third time step
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, 3);
    wd3 = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    sts_4th = EmeraldEarth.global_simulations!(gms, wd3, sts_3rd);

    # run the model for the fourth time step
    wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, 4);
    wd4 = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    sts_5th = EmeraldEarth.global_simulations!(gms, wd4, sts_4th);


    EmeraldEarth.global_simulations!(gms, wd4, sts_4th; single_thread = true, single_thread_regions = (:,111:180));
    begin
        i = 274;
        j = 112;
        gmij = gms[i,j];
        wdij = wd4[i,j];
        stij = deepcopy(sts_4th[i,j]);
        EmeraldEarth.grid_simulation!(gmij, wdij, stij);
    end;



    # use this commond to debug the code on a single thread and run the model grid by grid to find when the model crashes
    # then use the code block to debug that particular grid
    #=
    EmeraldEarth.global_simulations!(gms, wd2, sts_2nd; single_thread = true, single_thread_regions = (:,:));
    begin
        i = 305;
        j = 66;
        gmij = gms[i,j];
        wdij = wd2[i,j];
        stij = deepcopy(sts_2nd[i,j]);
        EmeraldEarth.grid_simulation!(gmij, wdij, stij);
    end;
    =#

    # run the global simulations for 10 time steps
    new_sts = sts_1st;
    for i in 1:10
        @time wds = EmeraldData.WeatherDrivers.preloaded_weather_drivers(WDRIVER_VER, 2019, 1, i);
        @time wdi = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
        @time new_sts = EmeraldEarth.global_simulations!(gms, wdi, new_sts);
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

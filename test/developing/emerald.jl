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
    wdi = EmeraldData.WeatherDrivers.weather_drivers_snapshot(wds, 1);
    mss = EmeraldData.WeatherDrivers.initial_states_snapshot(iss);

    # use multiple threads on curry and initialize the states matrix with multiple threads
    EmeraldEarth.add_threads!(32, FT);
    sts = EmeraldEarth.initial_states(gms, wdi, mss);

    # prepare the initial states matrix
    # TODO: create a function to load the first snapshot of the weather drivers (wdi) to save time




    EmeraldEarth.global_simulations!(gms, wdi, sts; single_thread = true);

    new_sts = EmeraldEarth.global_simulations!(gms, wdi, sts);



begin
    i = 1;
    j = 71;
    gmij = gms[i,j];
    wdij = wdi[i,j];
    ssij = mss[i,j];
    stij = deepcopy(sts[i,j]);
    EmeraldEarth.grid_simulation!(gmij, wdij, stij);
end;


    begin
        EmeraldLand.SPAC.initialize_spac!(EmeraldEarth.CACHE_CONFIG, EmeraldEarth.CACHE_SPAC, st11);
        EmeraldEarth.prescribe_traits_environment!(gm11, wd11);
        EmeraldLand.SPAC.spac!(EmeraldEarth.CACHE_CONFIG, EmeraldEarth.CACHE_SPAC, 3600);
    end;


    st11 = EmeraldEarth.grid_simulation!(gm11, wd11, st11);



    EmeraldEarth.grid_simulation!(gms[1,1], wdi[1,1], sts[1,1]);


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

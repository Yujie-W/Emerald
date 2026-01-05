# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using Profile

import Emerald.Land as ELAND
import GriddingMachine.Indexer as GMI

settings = ELAND.land_model_settings(mode = "testing");
settings["CONFIG_TAG"] = "testing_fvcb";
settings["MESSAGE_LEVEL"] = 0;
settings["REMOVE_WHEN_DONE"] = false;

gmd = GMI.grid_dict(GMI.LandDatasetLabels(settings["GM_VERSION"], 2019), 31.86389, 117.28083);
wd = GMI.grid_weather(GMI.WeatherDriverLabels("wd1",2019), 31.86389, 117.28083);
config = ELAND.site_config(settings);
spac = ELAND.site_spac(config, gmd);
sd = ELAND.parameters_to_save(settings["VARIABLES_TO_SAVE"]);
driver = ELAND.site_driver_tuple(gmd, wd);
results = ELAND.site_result_tuple(spac, wd, sd);

@time ELAND.simulation!(config, spac, driver, results, sd; selection = settings["SIMULATION_PERIOD"], δt = settings["TIME_STEP"]);

Profile.clear_malloc_data();

# 144 allocations due to a number of function calls to compute quantities such as GPP, ET, APAR, etc.
#   2 allocations due to the conversion of Dict{String,Any} to certain type such as Number
@time ELAND.simulation!(config, spac, driver, results, sd; selection = settings["SIMULATION_PERIOD"], δt = settings["TIME_STEP"]);


exit()

import Emerald.EmeraldFrontier as EF
import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldData.WeatherDrivers as WD
import Emerald.EmeraldLand.SPAC as SPAC

using Emerald.EmeraldIO.Jld2: read_jld2


# read in the data required for global simulations
# using Emerald.EmeraldIO.Jld2: read_jld2, save_jld2!
# jld2_to_read = "/mnt/net/ormosia/DATASERVER/model/Emerald/setups/emerald_grid_info_gm2_2019.jld2";
# jld_dicts = read_jld2(jld2_to_read, "GRID_INFO");
# for d in jld_dicts
#     if d["LAT_INDEX"] == 71 && d["LON_INDEX"] == 318
#         save_jld2!("debug.jld2", d);
#         break;
#     end
# end;


# load the data for debugging purpose
gmd = read_jld2("debug.jld2");
gmd["MESSAGE_LEVEL"] = 1;


# run the simulations with data from wd1
# df_result = EF.simulation!(wd_tag, gmd);
wd_tag = "wd1";

config = EF.spac_config(gmd);
config.ALLOW_LEAF_REGROWTH = false;
config.ALLOW_LEAF_SHEDDING = false;
config.ALLOW_XYLEM_GROWTH = false;
config.EFFECTIVE_LEAF_SPECTRA = false;
config.ENABLE_DROUGHT_LEGACY = false;
config.ENABLE_REF = true;
config.ENABLE_SIF = true;

spac = GD.grid_spac(config, gmd);
for s in spac.soils
    s.state.θ = s.trait.vc.Θ_SAT;
end;
spac.plant.pool.c_pool = Inf;
SPAC.initialize_spac!(config, spac);

df = WD.grid_weather_driver(wd_tag, gmd);
wdf = EF.prepare_wdf(spac, df);

EF.prescribe!(config, spac, wdf, 1; initialize_state = true);

spac_debug = deepcopy(spac);
idx_debug = 1;
for idx in eachindex(wdf.FDOY)
    @info "Debugging at time $idx" spac.plant.junction.state.v_storage spac.plant.junction.auxil.∂w∂t spac.plant.junction.auxil.pressure;
    global spac_debug = deepcopy(spac);
    global idx_debug = idx;
    try
        EF.simulation!(config, spac, wdf, idx);
    catch e
        @show e;
        break;
    end;
end;
spac_bak = deepcopy(spac_debug);
idx_bak = idx_debug;

#
#
# the step right before the error
# use the backed up spac to debug
#
#
spac = deepcopy(spac_bak);

for l in spac.plant.leaves
    l.flux.trait.g_limits[1] = 0;
end;

for idx in eachindex(wdf.FDOY)[idx_debug:end]
    @info "Debugging at time $idx" spac.plant.junction.state.v_storage spac.plant.junction.auxil.∂w∂t spac.plant.junction.auxil.pressure;
    global spac_debug = deepcopy(spac);
    global idx_debug = idx;
    try
        EF.simulation!(config, spac, wdf, idx);
    catch e
        @show e;
        break;
    end;
end;
spac_bak_2 = deepcopy(spac_debug);
idx_bak_2 = idx_debug;


# another round
spac = deepcopy(spac_bak_2); EF.simulation!(config, spac, wdf, idx_debug);

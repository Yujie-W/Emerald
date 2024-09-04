using Emerald
using Revise


FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);

EmeraldLand.SPAC.initialize_spac!(config, spac);
for i in 1:10
    @time EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @info "debugging" i spac.plant.pool.c_pool;
end;

spac.plant.pool.c_pool = -100;
for i in 1:10
    @time EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @info "debugging" i spac.plant.pool.c_pool;
end;

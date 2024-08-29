using Emerald
using Revise


FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);

EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(0));

for i in 1:500
    EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @info "debugging" i spac.plant.pool.c_pool;
end;

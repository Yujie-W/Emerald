using Emerald
using Revise


FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);

EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(0));
for i in 1:10
    EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    @info "debugging" i spac.plant.pool.c_pool;
end;

#
# Test the carbon pool use by new LAI growth
#
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1);
EmeraldLand.SPAC.spac!(config, spac, FT(3600));
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 5);
EmeraldLand.SPAC.spac!(config, spac, FT(3600));

#
# Test the comparison between the recovery and new growth schemes
#
spac.plant.branches[end].xylem.state.p_history .= -0.01;
EmeraldLand.SPAC.spac!(config, spac, FT(3600));
EmeraldLand.PlantHydraulics.xylem_recovery(spac.plant.branches[end].xylem, 10.0, spac.plant.branches[end].energy.s_aux.t);
spac.plant.branches[end].xylem.state.p_history .= -2;
EmeraldLand.SPAC.spac!(config, spac, FT(3600));
EmeraldLand.PlantHydraulics.xylem_recovery(spac.plant.branches[end].xylem, 10.0, spac.plant.branches[end].energy.s_aux.t);

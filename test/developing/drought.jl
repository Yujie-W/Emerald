using Emerald
using Revise


FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);

EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(1));


spac_bak1 = deepcopy(spac);
spac_bak2 = deepcopy(spac);
for i in 1:10
    global spac_bak1, spac_bak2;
    @time EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    spac_bak1 = deepcopy(spac);
    @time EmeraldLand.SPAC.spac!(config, spac, FT(3600));
    spac_bak2 = deepcopy(spac);
    # @info "debugging" spac_bak1.plant.leaves[1].xylem.auxil.pressure[end] spac_bak2.plant.leaves[1].xylem.auxil.pressure[end];
    # @info "debugging" [l.xylem.state.connected for l in spac.plant.leaves];
    @info "debugging" spac.plant.junction.state.v_storage spac.plant.junction.s_aux.pressure;
end;

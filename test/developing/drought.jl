using Emerald
using Revise


FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);

for s in spac.soils
    s.state.θ = s.trait.vc.Θ_RES + 0.001;
end;

EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(0));

while !spac.plant._leaf_shedded

    EmeraldLand.SPAC.spac!(config, spac, 3600);
    @info "Debugging" spac.plant.leaves[end].xylem.auxil.pressure[end] spac.plant.junction.state.v_storage spac.plant.junction.s_aux.pressure spac.plant.leaves[end].energy.s_aux.t;
    println();

end;

spac_bak = deepcopy(spac);




#  EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0);

begin "Debugging"
    spac = deepcopy(spac_bak);
    EmeraldLand.SPAC.prescribe_soil!(spac; swcs = (0.3, 0.3, 0.3, 0.3));
    EmeraldLand.SPAC.spac!(config, spac, 1);
    EmeraldLand.SPAC.spac!(config, spac, 3600);

    EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1);
    EmeraldLand.SPAC.spac!(config, spac, 3600);
end;




spac = deepcopy(spac_bak);
@info "Debugging" spac.plant.leaves[end].xylem.auxil.pressure[end] spac.plant.junction.state.v_storage spac.plant.junction.s_aux.pressure spac.plant.leaves[end].energy.s_aux.t;
EmeraldLand.SPAC.prescribe_soil!(spac; swcs = (0.3, 0.3, 0.3, 0.3));
EmeraldLand.SPAC.s_aux!(spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "Debugging" spac.plant.leaves[end].xylem.auxil.pressure[end] spac.plant.junction.state.v_storage spac.plant.junction.s_aux.pressure spac.plant.leaves[end].energy.s_aux.t;

EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1);
@info "Debugging" spac.plant.leaves[end].xylem.auxil.pressure[end] spac.plant.junction.state.v_storage spac.plant.junction.s_aux.pressure spac.plant.leaves[end].energy.s_aux.t;
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "Debugging" spac.plant.leaves[end].xylem.auxil.pressure[end] spac.plant.junction.state.v_storage spac.plant.junction.s_aux.pressure spac.plant.leaves[end].energy.s_aux.t;

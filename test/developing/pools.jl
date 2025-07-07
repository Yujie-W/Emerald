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




config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(1));
θ = 0.086;
for s in spac.soils
    s.state.θ = θ;
end;
psoil = EmeraldLand.SoilHydraulics.soil_ψ_25(spac.soils[1].trait.vc, θ) * EmeraldLand.PhysicalChemistry.relative_surface_tension(spac.soils[1].s_aux.t);
spac.plant.junction.state.v_storage = EmeraldLand.PlantHydraulics.capacitance_volume(spac.plant.junction.trait.pv, psoil, spac.plant.junction.s_aux.t) * spac.plant.junction.trait.v_max;
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@show θ EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.K_PLANT(spac) EmeraldLand.SPAC.K_PLANT(spac; include_leaf = false);

spac.plant.leaves[end].xylem.state.p_history
spac.plant.leaves[end].xylem.auxil.pressure

spac.plant.branches[end].xylem.state.p_history
spac.plant.branches[end].xylem.auxil.pressure

spac.plant.pool.c_pool = 10000;
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1, sai = 0, ci = 0.8);
spac.plant.leaves[end].xylem.state.p_history
spac.plant.leaves[end].xylem.auxil.pressure

EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 3, sai = 0, ci = 0.8);
spac.plant.leaves[end].xylem.state.p_history
spac.plant.leaves[end].xylem.auxil.pressure

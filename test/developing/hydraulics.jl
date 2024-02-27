import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.PlantHydraulics as PH


# steady state
config = NS.SPACConfiguration(Float64);
config.STEADY_STATE_FLOW = true;
spac = NS.BulkSPAC(config);
for leaf in spac.plant.leaves
    leaf.flux.state.g_H₂O_s_sunlit .= 0.1;
    leaf.g_H₂O_s_shaded = 0.1;
end;
PH.plant_water_budget!(spac, 1.0);
PH.plant_flow_profile!(config, spac);
PH.plant_pressure_profile!(spac);

for i in 1:1000
    PH.plant_water_budget!(spac, 2.0);
    PH.plant_flow_profile!(config, spac);
    PH.plant_pressure_profile!(spac);
    f_in = sum([PH.flow_out(root.xylem) for root in spac.plant.roots]);
    # @info "Flow out and into the junction" PH.flow_in(spac.plant.trunk.xylem) f_in spac.plant.junction.auxil.pressure spac.plant.junction.state.v_storage;
    @show spac.plant.junction.state.v_storage;
end;


# non-steady state
config.STEADY_STATE_FLOW = false;
spac = NS.BulkSPAC(config);
for leaf in spac.plant.leaves
    leaf.flux.state.g_H₂O_s_sunlit .= 0.1;
    leaf.g_H₂O_s_shaded = 0.1;
end;
PH.plant_water_budget!(spac, 1.0);
PH.plant_flow_profile!(config, spac);
PH.plant_pressure_profile!(spac);

for i in 1:1000
    PH.plant_water_budget!(spac, 2.0);
    PH.plant_flow_profile!(config, spac);
    PH.plant_pressure_profile!(spac);
    f_in = sum([PH.flow_out(root.xylem) for root in spac.plant.roots]);
    # @info "Flow out and into the junction" PH.flow_in(spac.plant.trunk.xylem) f_in spac.plant.junction.auxil.pressure spac.plant.junction.state.v_storage;
    @show spac.plant.junction.state.v_storage;
end;

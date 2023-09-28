import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.PlantHydraulics as PH


# steady state
config = NS.SPACConfiguration{Float64}(STEADY_STATE_FLOW = true);
spac = NS.MultiLayerSPAC(config);
for leaf in spac.LEAVES
    leaf.g_H₂O_s_sunlit .= 0.1;
    leaf.g_H₂O_s_shaded = 0.1;
end;
PH.plant_water_budget!(spac, 1.0);
PH.plant_flow_profile!(config, spac);
PH.plant_pressure_profile!(spac);

for i in 1:1000
    PH.plant_water_budget!(spac, 2.0);
    PH.plant_flow_profile!(config, spac);
    PH.plant_pressure_profile!(spac);
    f_in = sum([PH.flow_out(root.NS.xylem) for root in spac.ROOTS]);
    # @info "Flow out and into the junction" PH.flow_in(spac.TRUNK.NS.xylem) f_in spac.JUNCTION.auxil.pressure spac.JUNCTION.state.v_storage;
    @show spac.JUNCTION.state.v_storage;
end;


# non-steady state
config = NS.SPACConfiguration{Float64}(STEADY_STATE_FLOW = false);
spac = NS.MultiLayerSPAC(config);
for leaf in spac.LEAVES
    leaf.g_H₂O_s_sunlit .= 0.1;
    leaf.g_H₂O_s_shaded = 0.1;
end;
PH.plant_water_budget!(spac, 1.0);
PH.plant_flow_profile!(config, spac);
PH.plant_pressure_profile!(spac);

for i in 1:1000
    PH.plant_water_budget!(spac, 2.0);
    PH.plant_flow_profile!(config, spac);
    PH.plant_pressure_profile!(spac);
    f_in = sum([PH.flow_out(root.NS.xylem) for root in spac.ROOTS]);
    # @info "Flow out and into the junction" PH.flow_in(spac.TRUNK.NS.xylem) f_in spac.JUNCTION.auxil.pressure spac.JUNCTION.state.v_storage;
    @show spac.JUNCTION.state.v_storage;
end;

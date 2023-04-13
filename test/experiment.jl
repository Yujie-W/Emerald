#
# Experiment of whether to prescribe water and energy related parameters (already in runtests)
#
using Emerald;

wd_tag = "wd1";
gm_dict = EmeraldFrontier.gm_dict(EmeraldFrontier.GriddingMachineLabels(year=2019), 35.1, 115.2);
wd_data = EmeraldFrontier.weather_driver(wd_tag, gm_dict);
df_pres = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240, p_on = false, t_on = false, θ_on = false);
df_simu = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240, p_on = true, t_on = true, θ_on = true);


#
# Debugging the soil water budget
#
using Emerald;

FT = Float64;
config = EmeraldCore.Namespace.SPACConfiguration{FT}();
spac = EmeraldCore.Namespace.MultiLayerSPAC{FT}();
for slayer in spac.SOIL.LAYERS
    slayer.θ = 0.35;
end;
EmeraldCore.SPAC.initialize!(spac, config);
spac.METEO.rad_lw = 300;

spac.METEO.rain = 0;
for i in 1:10
    EmeraldCore.SPAC.spac!(spac, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_flow = sum([slayer.∂θ∂t * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_fout = spac.METEO.rain * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldCore.SoilHydraulics.root_sink(rlayer) for rlayer in spac.ROOTS]) * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldCore.SPAC.T_VEG(spac) * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    #@info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
end;

spac.METEO.rain = 1;
for i in 1:10
    EmeraldCore.SPAC.spac!(spac, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_flow = sum([slayer.∂θ∂t * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_fout = spac.METEO.rain * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldCore.SoilHydraulics.root_sink(rlayer) for rlayer in spac.ROOTS]) * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldCore.SPAC.T_VEG(spac) * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    #@info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
end;

spac.METEO.rad_sw.e_direct .= 0;
spac.METEO.rad_sw.e_diffuse .= 0;
spac.METEO.rain = 1;
for i in 1:10
    EmeraldCore.SPAC.spac!(spac, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_flow = sum([slayer.∂θ∂t * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_fout = spac.METEO.rain * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldCore.SoilHydraulics.root_sink(rlayer) for rlayer in spac.ROOTS]) * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldCore.SPAC.T_VEG(spac) * EmeraldCore.Constant.M_H₂O(FT) / EmeraldCore.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    @info "SWC per layer" [slayer.θ for slayer in spac.SOIL.LAYERS];
end;

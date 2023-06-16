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
config = EmeraldLand.Namespace.SPACConfiguration{FT}();
spac = EmeraldLand.Namespace.MultiLayerSPAC(config);
EmeraldLand.SPAC.update!(spac, config; swcs = (0.35, 0.35, 0.35, 0.35, 0.35));
EmeraldLand.SPAC.initialize!(spac, config);
spac.METEO.rad_lw = 300;

spac.METEO.rain = 0;
for i in 1:100
    EmeraldLand.SPAC.spac!(spac, config, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_flow = sum([slayer.∂θ∂t * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_fout = spac.METEO.rain * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldLand.SoilHydraulics.root_sink(rlayer) for rlayer in spac.ROOTS]) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldLand.SPAC.T_VEG(spac) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT);
    #@info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    #@info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
    pns = [slayer.TRACES.n_N₂ / (slayer.ΔZ * max(0, slayer.VC.Θ_SAT - slayer.θ)) / 1000 * slayer.t * EmeraldLand.Constant.GAS_R() for slayer in spac.SOIL.LAYERS];
    @info "N₂ partial pressure in kPa" pns;
    tss = [slayer.t for slayer in spac.SOIL.LAYERS];
    @info "Soil temperature in K" tss;
end;

spac.METEO.rain = 0.01;
for i in 1:10
    EmeraldLand.SPAC.spac!(spac, config, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_flow = sum([slayer.∂θ∂t * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_fout = spac.METEO.rain * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldLand.SoilHydraulics.root_sink(rlayer) for rlayer in spac.ROOTS]) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldLand.SPAC.T_VEG(spac) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    #@info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
end;

spac.METEO.rad_sw.e_direct .= 0;
spac.METEO.rad_sw.e_diffuse .= 0;
spac.METEO.rain = 0.01;
for i in 1:10
    EmeraldLand.SPAC.spac!(spac, config, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_flow = sum([slayer.∂θ∂t * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    soil_water_fout = spac.METEO.rain * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldLand.SoilHydraulics.root_sink(rlayer) for rlayer in spac.ROOTS]) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldLand.SPAC.T_VEG(spac) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    @info "SWC per layer" [slayer.θ for slayer in spac.SOIL.LAYERS];
end;




#
# Debugging LAI = 0, four steps
#     - RAD > 0    &&    LAI > 0
#     - RAD > 0    &&    LAI = 0
#     - RAD = 0    &&    LAI > 0
#     - RAD = 0    &&    LAI = 0
#
using Emerald;

function show_spac_info(node)
    beta = EmeraldLand.SPAC.BETA(spac);
    par = spac.CANOPY.RADIATION.par_in;
    ppar = EmeraldLand.SPAC.PPAR(spac);
    csif = EmeraldLand.SPAC.ΣSIF(spac);
    etr = EmeraldLand.SPAC.ΣETR(spac);
    gpp = EmeraldLand.SPAC.GPP(spac);
    ndvi = EmeraldLand.CanopyOptics.MODIS_NDVI(node);
    evi = EmeraldLand.CanopyOptics.MODIS_EVI(node);
    nirv = EmeraldLand.CanopyOptics.MODIS_NIRv(node);
    sifs = (EmeraldLand.CanopyOptics.TROPOMI_SIF683(node),
            EmeraldLand.CanopyOptics.TROPOMI_SIF740(node),
            EmeraldLand.CanopyOptics.OCO2_SIF759(node),
            EmeraldLand.CanopyOptics.OCO2_SIF770(node));
    tran = EmeraldLand.SPAC.T_VEG(spac);
    @info "SPAC Details" beta par ppar csif etr gpp tran ndvi evi nirv sifs;
end

FT = Float64;
config = EmeraldLand.Namespace.SPACConfiguration{FT}();
spac = EmeraldLand.Namespace.MultiLayerSPAC(config);
EmeraldLand.SPAC.update!(spac, config; swcs = (0.35, 0.35, 0.35, 0.35, 0.35));
EmeraldLand.SPAC.initialize!(spac, config);
spac.METEO.rad_lw = 300;

@info "RAD > 0 and LAI = 0";
EmeraldLand.SPAC.update!(spac, config; lai = 0);
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD > 0 and LAI > 0";
EmeraldLand.SPAC.update!(spac, config; lai = 1);
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD > 0 and LAI = 0";
EmeraldLand.SPAC.update!(spac, config; lai = 0);
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD > 0 and LAI > 0";
EmeraldLand.SPAC.update!(spac, config; lai = 1);
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI = 0";
EmeraldLand.SPAC.update!(spac, config; lai = 0);
spac.METEO.rad_sw.e_direct .= 0;
spac.METEO.rad_sw.e_diffuse .= 0;
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI > 0";
EmeraldLand.SPAC.update!(spac, config; lai = 1);
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI > 0 (SZA > 90)";
spac.ANGLES.sza = 90;
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI = 0 (SZA > 90)";
EmeraldLand.SPAC.update!(spac, config; lai = 0);
EmeraldLand.SPAC.spac!(spac, config, FT(360));
show_spac_info(spac);

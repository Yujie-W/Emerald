#
# Experiment of whether to prescribe water and energy related parameters (already in runtests)
#
using Emerald;

wd_tag = "wd1";
gm_dict = EmeraldFrontier.gm_dict(EmeraldFrontier.GriddingMachineLabels(year=2019), 35.1, 115.2);
wd_data = EmeraldFrontier.weather_driver(wd_tag, gm_dict);
df_pres = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240);
df_simu = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240);


#
# Debugging the root disconnection and reconnection issues
#
using Emerald;

FT = Float64;
config = EmeraldLand.Namespace.SPACConfiguration{FT}(DEBUG = true);
spac = EmeraldLand.Namespace.BulkSPAC(config);
spac.meteo.rad_lw = 300;
EmeraldLand.SPAC.initialize!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(360));


#
# Debugging the soil water budget
#
using Emerald;

FT = Float64;
config = EmeraldLand.Namespace.SPACConfiguration{FT}(DEBUG = true);
spac = EmeraldLand.Namespace.BulkSPAC(config);
#EmeraldLand.SPAC.prescribe_soil!(spac; swcs = (0.35, 0.35, 0.43, 0.35, 0.43));
EmeraldLand.SPAC.initialize!(config, spac);
spac.meteo.rad_lw = 300;
EmeraldLand.SPAC.spac!(config, spac, FT(360));

spac.meteo.rain = 0;
for i in 1:100
    EmeraldLand.SPAC.spac!(config, spac, FT(360));
    tswc = sum([soil.θ * soil.auxil.δz for soil in spac.SOIL.LAYERS]);
    soil_water_flow = sum([soil.∂θ∂t * soil.auxil.δz for soil in spac.SOIL.LAYERS]);
    soil_water_fout = spac.meteo.rain * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldLand.SoilHydraulics.root_sink(rlayer) for rlayer in spac.plant.roots]) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldLand.SPAC.T_VEG(spac) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT);
    #@info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    #@info "Total water is" tswc [soil.θ for soil in spac.SOIL.LAYERS];
    pns = [soil.state.ns[4] / (soil.auxil.δz * max(0, soil.VC.Θ_SAT - soil.θ)) / 1000 * soil.t * EmeraldLand.Constant.GAS_R() for soil in spac.SOIL.LAYERS];
    @info "N₂ partial pressure in kPa" pns;
    tss = [soil.t for soil in spac.SOIL.LAYERS];
    @info "Soil temperature in K" tss;
end;

spac.meteo.rain = 0.01;
for i in 1:10
    EmeraldLand.SPAC.spac!(config, spac, FT(360));
    tswc = sum([soil.θ * soil.auxil.δz for soil in spac.SOIL.LAYERS]);
    soil_water_flow = sum([soil.∂θ∂t * soil.auxil.δz for soil in spac.SOIL.LAYERS]);
    soil_water_fout = spac.meteo.rain * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldLand.SoilHydraulics.root_sink(rlayer) for rlayer in spac.plant.roots]) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldLand.SPAC.T_VEG(spac) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    #@info "Total water is" tswc [soil.θ for soil in spac.SOIL.LAYERS];
end;

spac.meteo.rad_sw.e_dir .= 0;
spac.meteo.rad_sw.e_dif .= 0;
spac.meteo.rain = 0.01;
for i in 1:10
    EmeraldLand.SPAC.spac!(config, spac, FT(360));
    tswc = sum([soil.θ * soil.auxil.δz for soil in spac.SOIL.LAYERS]);
    soil_water_flow = sum([soil.∂θ∂t * soil.auxil.δz for soil in spac.SOIL.LAYERS]);
    soil_water_fout = spac.meteo.rain * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) - soil_water_flow;
    root_water_flow = sum([EmeraldLand.SoilHydraulics.root_sink(rlayer) for rlayer in spac.plant.roots]) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT) / spac.SOIL.AREA;
    leaf_water_flow = EmeraldLand.SPAC.T_VEG(spac) * EmeraldLand.Constant.M_H₂O(FT) / EmeraldLand.Constant.ρ_H₂O(FT);
    @info "Debugging" soil_water_flow soil_water_fout root_water_flow leaf_water_flow spac.SOIL.runoff;
    @info "SWC per layer" [soil.θ for soil in spac.SOIL.LAYERS];
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
    par = EmeraldLand.SPAC.PAR(config, spac);
    ppar = EmeraldLand.SPAC.PPAR(spac);
    csif = EmeraldLand.SPAC.ΣSIF(spac);
    etr = EmeraldLand.SPAC.ΣETR(spac);
    gpp = EmeraldLand.SPAC.GPP(spac);
    ndvi = EmeraldLand.CanopyOptics.MODIS_NDVI(node);
    evi = EmeraldLand.CanopyOptics.MODIS_EVI(node);
    nirv = EmeraldLand.CanopyOptics.MODIS_NIRv(node);
    sifs = (EmeraldLand.CanopyOptics.TROPOMI_SIF683(node), EmeraldLand.CanopyOptics.TROPOMI_SIF740(node), EmeraldLand.CanopyOptics.OCO2_SIF759(node), EmeraldLand.CanopyOptics.OCO2_SIF770(node));
    tran = EmeraldLand.SPAC.T_VEG(spac);
    @info "SPAC Details" beta par ppar csif etr gpp tran ndvi evi nirv sifs;
end;

FT = Float64;
config = EmeraldLand.Namespace.SPACConfiguration{FT}();
spac = EmeraldLand.Namespace.BulkSPAC(config);
EmeraldLand.SPAC.prescribe_soil!(spac; swcs = (0.35, 0.35, 0.35, 0.35, 0.35));
EmeraldLand.SPAC.initialize!(config, spac);
spac.meteo.rad_lw = 300;

@info "RAD > 0 and LAI = 0";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0);
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD > 0 and LAI > 0";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1);
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD > 0 and LAI = 0";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0);
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD > 0 and LAI > 0";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1);
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI = 0";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0);
spac.meteo.rad_sw.e_dir .= 0;
spac.meteo.rad_sw.e_dif .= 0;
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI > 0";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 1);
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI > 0 (SZA > 90)";
spac.canopy.sun_geometry.state.sza = 90;
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

@info "RAD = 0 and LAI = 0 (SZA > 90)";
EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0);
EmeraldLand.SPAC.spac!(config, spac, FT(360));
show_spac_info(spac);

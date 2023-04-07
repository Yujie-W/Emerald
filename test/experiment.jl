#
# Experiment of whether to prescribe water and energy related parameters
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
spac = EmeraldCore.Namespace.MultiLayerSPAC{FT}();
for slayer in spac.SOIL.LAYERS
    slayer.θ = 0.35;
end;
EmeraldCore.SPAC.initialize!(spac);
spac.RAD_LW = 300;

spac.METEO.rain = 0;
for i in 1:10
    EmeraldCore.SPAC.spac!(spac, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    @info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
end;

spac.METEO.rain = 10;
for i in 1:10
    EmeraldCore.SPAC.spac!(spac, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    @info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
    @info "Runoff" spac.SOIL.runoff;
end;

spac.RAD_SW.e_direct .= 0;
spac.RAD_SW.e_diffuse .= 0;
spac.METEO.rain = 10;
for i in 1:10
    EmeraldCore.SPAC.spac!(spac, FT(360));
    tswc = sum([slayer.θ * slayer.ΔZ for slayer in spac.SOIL.LAYERS]);
    @info "Total water is" tswc [slayer.θ for slayer in spac.SOIL.LAYERS];
    @info "Runoff" spac.SOIL.runoff;
end;

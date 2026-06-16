using DataFrames
using Emerald
using Emerald.Namespace
using Emerald.SPAC
using PkgUtility.UniversalConstants: energy_to_photon


FT = Float64;

config = Emerald.Namespace.SPACConfig(FT);
#config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTSCOPE();
config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTEmerald();

lai = 3;

spac = Namespace.BulkSPAC(config; air_bounds = collect(0:0.25:13));

#spac.canopy.structure.trait.lidf.A = 0;
#spac.canopy.structure.trait.lidf.B = 1;

# LEAF:CHL SIF escape is WL-dependent because of the SIF spectrum shape
wl = 500;
mask = (wl - 10) .<= config.CONSTANTS.SPECTRA.Λ .<= (wl + 10);
spac.meteo.rad_sw.e_dir[.!mask] .= 0;
spac.meteo.rad_sw.e_dif[.!mask] .= 0;
spac.meteo.rad_sw.e_dir .= 0;
#spac.meteo.rad_sw.e_dif .= 0;

spac.canopy.sun_geometry.state.sza = 30;
spac.canopy.sun_geometry.state.saa = 180;
spac.canopy.sensor_geometry.state.vza = 0;
spac.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(config, spac; sai = 0, lai = lai);
SPAC.initialize_spac!(config, spac);
SPAC.spac!(config, spac, 1);


f_sife = spac.plant.leaves[1].bio.auxil.f_sife[config.CONSTANTS.SPECTRA.IΛ_SIFE];
e_net = (spac.canopy.sun_geometry.auxil.e_net_dif .+ spac.canopy.sun_geometry.auxil.e_net_dir)[config.CONSTANTS.SPECTRA.IΛ_SIFE, :] .* f_sife;
sif_chl = deepcopy(spac.canopy.sun_geometry.auxil.e_sif_chl);
sif_leaf = spac.canopy.sun_geometry.auxil.e_sifꜛ_layer .+ spac.canopy.sun_geometry.auxil.e_sifꜜ_layer;

e_net_sum = [config.CONSTANTS.SPECTRA.ΔΛ_SIFE' * energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIFE,eee) * 1e6 for eee in eachcol(e_net)];
sif_chl_sum = [config.CONSTANTS.SPECTRA.ΔΛ_SIF' * energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIF,sss) * 1e6 for sss in eachcol(sif_chl)];
sif_leaf_sum = [config.CONSTANTS.SPECTRA.ΔΛ_SIF' * energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIF,sss) * 1e6 for sss in eachcol(sif_leaf)];

@show sif_chl_sum;
@show e_net_sum;

phi_layer = sif_chl_sum ./ e_net_sum;
esc_leaf = sif_leaf_sum ./ sif_chl_sum;
for i in eachindex(phi_layer)
    @show i,phi_layer[i],esc_leaf[i];
end;


# double check if total SW radiation is conserved
sw_in = config.CONSTANTS.SPECTRA.ΔΛ' * (spac.meteo.rad_sw.e_dir .+ spac.meteo.rad_sw.e_dif) / 1000
sw_out = config.CONSTANTS.SPECTRA.ΔΛ' * spac.canopy.sun_geometry.auxil.e_difꜛ[:,1] / 1000
net_sw_plant = sum(config.CONSTANTS.SPECTRA.ΔΛ' * (spac.canopy.sun_geometry.auxil.e_net_dif .+ spac.canopy.sun_geometry.auxil.e_net_dir) ./ 1000)
net_sw_soil = sum(config.CONSTANTS.SPECTRA.ΔΛ' * (spac.soil_bulk.auxil.e_net_dif .+ spac.soil_bulk.auxil.e_net_dir) ./ 1000)

@info "Energy balance check" sw_in sw_out + net_sw_plant + net_sw_soil;

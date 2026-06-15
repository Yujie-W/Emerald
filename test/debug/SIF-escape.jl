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

spac.canopy.sun_geometry.state.sza = 30;
spac.canopy.sun_geometry.state.saa = 180;
spac.canopy.sensor_geometry.state.vza = 0;
spac.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(config, spac; sai = 0, lai = lai);
SPAC.initialize_spac!(config, spac);
SPAC.spac!(config, spac, 1);


e_net = (spac.canopy.sun_geometry.auxil.e_net_dif .+ spac.canopy.sun_geometry.auxil.e_net_dir)[config.CONSTANTS.SPECTRA.IΛ_SIFE, :];
sif_chl = deepcopy(spac.canopy.sun_geometry.auxil.e_sif_chl);
sif_leaf = spac.canopy.sun_geometry.auxil.e_sifꜛ_layer .+ spac.canopy.sun_geometry.auxil.e_sifꜜ_layer;

e_net_sum = [config.CONSTANTS.SPECTRA.ΔΛ_SIFE' * energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIFE,eee) * 1e6 for eee in eachcol(e_net)];
sif_chl_sum = [config.CONSTANTS.SPECTRA.ΔΛ_SIF' * energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIF,sss) * 1e6 for sss in eachcol(sif_chl)];
sif_leaf_sum = [config.CONSTANTS.SPECTRA.ΔΛ_SIF' * energy_to_photon.(config.CONSTANTS.SPECTRA.Λ_SIF,sss) * 1e6 for sss in eachcol(sif_leaf)];

phi_layer = sif_chl_sum ./ e_net_sum;
esc_leaf = sif_leaf_sum ./ sif_chl_sum;
for i in eachindex(phi_layer)
    @show i,phi_layer[i],esc_leaf[i];
end;

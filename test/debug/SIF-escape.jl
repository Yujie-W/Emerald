using DataFrames
using Emerald
using Emerald.Namespace
using Emerald.SPAC


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


sif_chl = spac.canopy.sun_geometry.auxil.e_sif_chl;
sif_leaf = spac.canopy.sun_geometry.auxil.e_sifꜛ_layer .+ spac.canopy.sun_geometry.auxil.e_sifꜜ_layer;

esc_leaf = sif_leaf ./ sif_chl

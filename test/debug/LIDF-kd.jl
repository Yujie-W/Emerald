using DataFrames
using Emerald
using Emerald.Namespace
using Emerald.SPAC


FT = Float64;

config = Emerald.Namespace.SPACConfig(FT);
#config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTEmerald();
config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTSCOPE();
lai = 3;
spac = Namespace.BulkSPAC(config; air_bounds = collect(0:0.25:13));
spac.canopy.sun_geometry.state.sza = 30;
spac.canopy.sun_geometry.state.saa = 180;
spac.canopy.sensor_geometry.state.vza = 0;
spac.canopy.sensor_geometry.state.vaa = 0;
spac.meteo.rad_sw.e_dir .= 0;
SPAC.prescribe_traits!(config, spac; sai = 0, lai = lai);

SPAC.initialize_spac!(config, spac);
SPAC.spac!(config, spac, 1);

# 1, 5, 6
lidfs = [[0,0], [-0.35,-0.15], [-1,0], [1,0], [0,-1], [0,1]];
for b in -1:0.1:1
    spac.canopy.structure.trait.lidf.A = 0;
    spac.canopy.structure.trait.lidf.B = b;
    SPAC.t_aux!(config, spac);
end;
for a in -1:0.1:1
    spac.canopy.structure.trait.lidf.A = a;
    spac.canopy.structure.trait.lidf.B = 0;
    SPAC.t_aux!(config, spac);
end;
for l in lidfs
    spac.canopy.structure.trait.lidf.A = l[1];
    spac.canopy.structure.trait.lidf.B = l[2];
    SPAC.t_aux!(config, spac);
end;


#spac.canopy.structure.trait.lidf.A = 0;
#spac.canopy.structure.trait.lidf.B = 1;

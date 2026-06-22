using DataFrames
using Emerald
using Emerald.Namespace
using Emerald.SPAC


FT = Float64;


#=
config = Emerald.Namespace.SPACConfig(FT);
config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTEmerald();
#config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTSCOPE();
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
    spac.canopy.structure.trait.lidf.A = b;
    spac.canopy.structure.trait.lidf.B = 0;
    SPAC.t_aux!(config, spac);
    SPAC.spac!(config, spac, 0);
    @show b;
    # @show spac.canopy.structure.auxil.p_incl_leaf;
    @show spac.canopy.structure.auxil.τ_dd_diffuse[1];
    @show spac.canopy.structure.auxil.bf_leaf;
    @show spac.canopy.sun_geometry.auxil.e_difꜛ[6:10,1];
end;
=#

#spac.canopy.structure.trait.lidf.A = 0;
#spac.canopy.structure.trait.lidf.B = 1;


config = Emerald.Namespace.SPACConfig(FT);
#config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTEmerald();
config.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTSCOPE();
lai = 3;

for nl in 5:200
    dair = 6 / nl;
    air_bounds = collect(0:dair:13);
    spac = Namespace.BulkSPAC(config; air_bounds);
    spac.canopy.sun_geometry.state.sza = 0;
    spac.canopy.sun_geometry.state.saa = 180;
    spac.canopy.sensor_geometry.state.vza = 0;
    spac.canopy.sensor_geometry.state.vaa = 0;
    spac.meteo.rad_sw.e_dif .= 0;
    SPAC.prescribe_traits!(config, spac; sai = 0, lai = lai);
    SPAC.initialize_spac!(config, spac);
    SPAC.spac!(config, spac, 1);
    println(spac.canopy.sun_geometry.auxil.e_difꜛ[26,1]);
end;

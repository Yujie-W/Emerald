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
spac.canopy.sun_geometry.state.sza = 30;
spac.canopy.sun_geometry.state.saa = 180;
spac.canopy.sensor_geometry.state.vza = 0;
spac.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(config, spac; sai = 0, lai = lai);
SPAC.initialize_spac!(config, spac);

# SPAC.spac!(config, spac, 1);


@info "geometry";
for angle in collect(FT, -89.5:0.5:89.5)
    spac.canopy.sensor_geometry.state.vza = abs(angle);
    if angle <= 0
        spac.canopy.sensor_geometry.state.vaa = 0;
    else
        spac.canopy.sensor_geometry.state.vaa = 180;
    end;

    SPAC.spac!(config, spac, 0);

    println(angle, ",", spac.canopy.sensor_geometry.auxil.reflectance[54], ",", spac.canopy.sun_geometry.auxil.albedo[54]);
end;


#=
spac.canopy.sensor_geometry.state.vza = 74;
spac.canopy.sensor_geometry.state.vaa = 0;
SPAC.spac!(config, spac, 0);
println(74, ",", spac.canopy.sensor_geometry.auxil.reflectance[54]);

spac.canopy.sensor_geometry.state.vza = 65;
spac.canopy.sensor_geometry.state.vaa = 0;
SPAC.spac!(config, spac, 0);
println(65, ",", spac.canopy.sensor_geometry.auxil.reflectance[54]);
=#

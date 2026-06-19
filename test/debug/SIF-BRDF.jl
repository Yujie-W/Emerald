using DataFrames
using Emerald
using Emerald.Namespace
using Emerald.SPAC


FT = Float64;

configs = Emerald.Namespace.SPACConfig(FT);
configs.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTSCOPE();
configs.METHODS.SOIL_ALBEDO = Namespace.SoilAlbedoPrescribe();
confige = Emerald.Namespace.SPACConfig(FT);
confige.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTEmerald();
confige.METHODS.SOIL_ALBEDO = Namespace.SoilAlbedoPrescribe();

lai = 3;

spacs = Namespace.BulkSPAC(configs; air_bounds = collect(0:0.25:13));
spacs.soil_bulk.auxil.ρ_sw .= 0.2;
spacs.canopy.sun_geometry.state.sza = 30;
spacs.canopy.sun_geometry.state.saa = 180;
spacs.canopy.sensor_geometry.state.vza = 0;
spacs.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(configs, spacs; sai = 0, lai = lai);
SPAC.initialize_spac!(configs, spacs);

space = Namespace.BulkSPAC(confige; air_bounds = collect(0:0.25:13));
space.soil_bulk.auxil.ρ_sw .= 0.2;
space.canopy.sun_geometry.state.sza = 30;
space.canopy.sun_geometry.state.saa = 180;
space.canopy.sensor_geometry.state.vza = 0;
space.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(confige, space; sai = 0, lai = lai);
SPAC.initialize_spac!(confige, space);


SPAC.spac!(configs, spacs, 1);
SPAC.spac!(confige, space, 1);


@info "SCOPE and Emerald SIF geometry";
for angle in collect(FT, -89.5:0.5:89.5)
    spacs.canopy.sensor_geometry.state.vza = abs(angle);
    space.canopy.sensor_geometry.state.vza = abs(angle);
    if angle <= 0
        spacs.canopy.sensor_geometry.state.vaa = 0;
        space.canopy.sensor_geometry.state.vaa = 0;
    else
        spacs.canopy.sensor_geometry.state.vaa = 180;
        space.canopy.sensor_geometry.state.vaa = 180;
    end;

    SPAC.spac!(configs, spacs, 0);
    SPAC.spac!(confige, space, 0);

    println(angle, ",", spacs.canopy.sensor_geometry.auxil.sif_obs[20], ",", space.canopy.sensor_geometry.auxil.sif_obs[20]);
end;

using DataFrames
using Emerald
using Emerald.Namespace
using Emerald.SPAC


FT = Float64;

config_scope = Emerald.Namespace.SPACConfig(FT);
config_scope.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTSCOPE();
config_emerald = Emerald.Namespace.SPACConfig(FT);
config_emerald.METHODS.CANOPY_RT_METHOD = Namespace.CanopyRTEmerald();

lai = 3;

spac_scope = Namespace.BulkSPAC(config_scope);
spac_scope.canopy.sun_geometry.state.sza = 30;
spac_scope.canopy.sun_geometry.state.saa = 180;
spac_scope.canopy.sensor_geometry.state.vza = 80;
spac_scope.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(config_scope, spac_scope; sai = 0, lai = lai);
SPAC.initialize_spac!(config_scope, spac_scope);

spac_emerald = Namespace.BulkSPAC(config_emerald);
spac_emerald.canopy.sun_geometry.state.sza = 30;
spac_emerald.canopy.sun_geometry.state.saa = 180;
spac_emerald.canopy.sensor_geometry.state.vza = 80;
spac_emerald.canopy.sensor_geometry.state.vaa = 0;
SPAC.prescribe_traits!(config_emerald, spac_emerald; sai = 0, lai = lai);
SPAC.initialize_spac!(config_emerald, spac_emerald);

SPAC.spac!(config_scope, spac_scope, 1);
SPAC.spac!(config_emerald, spac_emerald, 1);

#@info "debugging SCOPE" spac_scope.canopy.structure.auxil.w_ddb_leaf spac_scope.canopy.structure.auxil.w_ddf_leaf spac_scope.canopy.structure.auxil.w_ddb_stem spac_scope.canopy.structure.auxil.w_ddf_stem spac_scope.canopy.sun_geometry.auxil.w_sdb_leaf spac_scope.canopy.sun_geometry.auxil.w_sdf_leaf spac_scope.canopy.sun_geometry.auxil.w_sdb_stem spac_scope.canopy.sun_geometry.auxil.w_sdf_stem spac_scope.canopy.sensor_geometry.auxil.w_dob_leaf spac_scope.canopy.sensor_geometry.auxil.w_dof_leaf spac_scope.canopy.sensor_geometry.auxil.w_dob_stem spac_scope.canopy.sensor_geometry.auxil.w_dof_stem;
#@info "debugging Emerald" spac_emerald.canopy.structure.auxil.w_ddb_leaf spac_emerald.canopy.structure.auxil.w_ddf_leaf spac_emerald.canopy.structure.auxil.w_ddb_stem spac_emerald.canopy.structure.auxil.w_ddf_stem spac_emerald.canopy.sun_geometry.auxil.w_sdb_leaf spac_emerald.canopy.sun_geometry.auxil.w_sdf_leaf spac_emerald.canopy.sun_geometry.auxil.w_sdb_stem spac_emerald.canopy.sun_geometry.auxil.w_sdf_stem spac_emerald.canopy.sensor_geometry.auxil.w_dob_leaf spac_emerald.canopy.sensor_geometry.auxil.w_dof_leaf spac_emerald.canopy.sensor_geometry.auxil.w_dob_stem spac_emerald.canopy.sensor_geometry.auxil.w_dof_stem;



@info "spectra";
r_scope = spac_scope.canopy.sun_geometry.auxil.e_difꜛ[:,1] ./ (spac_scope.meteo.rad_sw.e_dif .+ spac_scope.meteo.rad_sw.e_dir);
r_emerald = spac_emerald.canopy.sun_geometry.auxil.e_difꜛ[:,1] ./ (spac_emerald.meteo.rad_sw.e_dif .+ spac_emerald.meteo.rad_sw.e_dir);
brf_scope = spac_scope.canopy.sensor_geometry.auxil.reflectance;
brf_emerald = spac_emerald.canopy.sensor_geometry.auxil.reflectance;
for i in eachindex(spac_emerald.canopy.sensor_geometry.auxil.reflectance)
    println(config_scope.CONSTANTS.SPECTRA.Λ[i], ",", r_scope[i], ",", r_emerald[i], ",", brf_scope[i], ",", brf_emerald[i]);
end;


#=
@info "geometry";
for angle in collect(-87.5:5:90)
    spac_scope.canopy.sensor_geometry.state.vza = abs(angle);
    spac_emerald.canopy.sensor_geometry.state.vza = abs(angle);
    if angle < 0
        spac_scope.canopy.sensor_geometry.state.vaa = 0;
        spac_emerald.canopy.sensor_geometry.state.vaa = 0;
    else
        spac_scope.canopy.sensor_geometry.state.vaa = 180;
        spac_emerald.canopy.sensor_geometry.state.vaa = 180;
    end;

    SPAC.spac!(config_scope, spac_scope, 1);
    SPAC.spac!(config_emerald, spac_emerald, 1);

    @info "angle" angle;
    #@info "debugging SCOPE" spac_scope.canopy.structure.auxil.w_ddb_leaf spac_scope.canopy.structure.auxil.w_ddf_leaf spac_scope.canopy.structure.auxil.w_ddb_stem spac_scope.canopy.structure.auxil.w_ddf_stem spac_scope.canopy.sun_geometry.auxil.w_sdb_leaf spac_scope.canopy.sun_geometry.auxil.w_sdf_leaf spac_scope.canopy.sun_geometry.auxil.w_sdb_stem spac_scope.canopy.sun_geometry.auxil.w_sdf_stem spac_scope.canopy.sensor_geometry.auxil.w_dob_leaf spac_scope.canopy.sensor_geometry.auxil.w_dof_leaf spac_scope.canopy.sensor_geometry.auxil.w_dob_stem spac_scope.canopy.sensor_geometry.auxil.w_dof_stem;
    #@info "debugging Emerald" spac_emerald.canopy.structure.auxil.w_ddb_leaf spac_emerald.canopy.structure.auxil.w_ddf_leaf spac_emerald.canopy.structure.auxil.w_ddb_stem spac_emerald.canopy.structure.auxil.w_ddf_stem spac_emerald.canopy.sun_geometry.auxil.w_sdb_leaf spac_emerald.canopy.sun_geometry.auxil.w_sdf_leaf spac_emerald.canopy.sun_geometry.auxil.w_sdb_stem spac_emerald.canopy.sun_geometry.auxil.w_sdf_stem spac_emerald.canopy.sensor_geometry.auxil.w_dob_leaf spac_emerald.canopy.sensor_geometry.auxil.w_dof_leaf spac_emerald.canopy.sensor_geometry.auxil.w_dob_stem spac_emerald.canopy.sensor_geometry.auxil.w_dof_stem;

    #println(angle, ",", spac_scope.canopy.sensor_geometry.auxil.reflectance[54], ",", spac_emerald.canopy.sensor_geometry.auxil.reflectance[54]);
end;
=#

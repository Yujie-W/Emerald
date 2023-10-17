# This file contains functions to compute the canopy optical properties at the sensor direction

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-12: add function sensor_spectrum! to compute the spectra at the sensor
#     2023-Oct-14: do nothing if REF is not enabled
#     2023-Oct-14: if SZA > 89, set all shortwave fluxes to 0 and reflectance to NaN
#     2023-Oct-14: if LAI <= 0, use soil reflectance only
#
#######################################################################################################################################################################################################
"""

    reflection_spectrum!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Computes the spectra at the sensor direction, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function reflection_spectrum!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    if !config.ENABLE_REF
        return nothing
    end;

    (; DIM_LAYER, DIM_WL) = config;
    can_str = spac.canopy.structure;
    rad_sw = spac.meteo.rad_sw;
    sbulk = spac.soil_bulk;
    sen_geo = spac.canopy.sensor_geometry;
    sun_geo = spac.canopy.sun_geometry;

    if sun_geo.state.sza > 89
        sen_geo.auxil.e_sensor_layer .= 0;
        sen_geo.auxil.e_sensor .= 0;
        sen_geo.auxil.reflectance .= NaN;

        return nothing
    end;

    if can_str.state.lai <= 0
        sen_geo.auxil.e_sensor_layer .= 0;
        sen_geo.auxil.e_sensor_layer[:,end] .= view(sun_geo.auxil.e_difꜛ,:,DIM_LAYER+1);
        sen_geo.auxil.e_sensor .= view(sun_geo.auxil.e_difꜛ,:,DIM_LAYER+1) ./ FT(π);
        sen_geo.auxil.reflectance .= sbulk.auxil.ρ_sw;

        return nothing
    end;

    # Run the canopy optical properties simulations only if canopy reflectance feature is enabled

    # compute the spectra at the observer direction
    for i in 1:DIM_LAYER
        e_d_i = view(sun_geo.auxil.e_difꜜ,:,i);             # downward diffuse radiation at upper boundary
        e_u_i = view(sun_geo.auxil.e_difꜛ,:,i);             # upward diffuse radiation at upper boundary
        sen_i = view(sen_geo.auxil.e_sensor_layer,:,i);  # radiation towards the viewing direction per layer (including soil)

        dob_i = view(sen_geo.auxil.dob_leaf,:,i);        # scattering coefficient backward for diffuse->observer
        dof_i = view(sen_geo.auxil.dof_leaf,:,i);        # scattering coefficient forward for diffuse->observer
        so_i  = view(sen_geo.auxil.so_leaf ,:,i);        # bidirectional from solar to observer

        sen_i  .= sen_geo.auxil.p_sensor[i]     .* dob_i .* e_d_i .+
                  sen_geo.auxil.p_sensor[i]     .* dof_i .* e_u_i .+
                  sen_geo.auxil.p_sun_sensor[i] .* so_i  .* rad_sw.e_dir;
        sen_i .*= can_str.state.δlai[i] * can_str.auxil.ci;
    end;
    sen_geo.auxil.e_sensor_layer[:,end] .= sen_geo.auxil.p_sensor_soil .* view(sun_geo.auxil.e_difꜛ,:,DIM_LAYER+1);

    # compute the spectra at the sensor
    for i in 1:DIM_WL
        sen_geo.auxil.e_sensor[i] = sum(view(sen_geo.auxil.e_sensor_layer,i,:)) / FT(π);
    end;

    # Note, this reflectance calculation is not correct because the sun-sensor geometry is not taken into account (reflectance is not isotropic)
    #     This is to compare with remote sensing data, which use the same calculation (* π)
    #     SCOPE does this a bit differently (controlling numerical issues), but still not correct
    # For real albedo, which is the surface reflectance, it needs to be the ratio between upward diffuse light and total downward light
    sen_geo.auxil.reflectance .= sen_geo.auxil.e_sensor .* FT(π) ./ (rad_sw.e_dir .+ rad_sw.e_dif);

    return nothing
end;

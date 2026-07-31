# This file contains functions to compute the canopy optical properties at the sensor direction

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-12: add function sensor_spectrum! to compute the spectra at the sensor
#     2023-Oct-14: do nothing if REF is not enabled
#     2023-Oct-14: if SZA > 89, set all shortwave fluxes to 0 and reflectance to NaN
#     2023-Oct-14: if LAI <= 0, use soil reflectance only
#     2023-Oct-18: account for SAI in the canopy reflectance calculation
#     2024-Sep-04: separate leaf and stem optical properties
# Bug fixes
#     2024-Mar-06: ci impact on fraction from viewer direction (otherwise will be accounted twice)
#     2025-Sep-12: add a special case when toral rad is zero (to avoid NaN issue)
#
#######################################################################################################################################################################################################
"""

    reflection_spectrum!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Computes the spectra at the sensor direction, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function reflection_spectrum!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    if !config.FEATURES.ENABLE_REF
        return nothing
    end;

    can_str = spac.canopy.structure;
    rad_sw = spac.meteo.rad_sw;
    sbulk = spac.soil_bulk;
    sen_geo = spac.canopy.sensor_geometry;
    sun_geo = spac.canopy.sun_geometry;
    n_layer = length(can_str.trait.δlai);
    rad_sw = spac.meteo.rad_sw;
    (; SPECTRA) = config.CONSTANTS;

    # if sza > 89, set all the radiation variables to 0
    total_sw_rad = (rad_sw.e_dir' * SPECTRA.ΔΛ + rad_sw.e_dif' * SPECTRA.ΔΛ) / 1000;
    if sun_geo.state.sza > 89 || total_sw_rad <= 0
        sen_geo.auxil.e_sensor_layer .= 0;
        sen_geo.auxil.e_sensor .= 0;
        sen_geo.auxil.reflectance .= NaN;

        return nothing
    end;

    if can_str.trait.lai <= 0 && can_str.trait.sai <= 0
        sen_geo.auxil.e_sensor_layer .= 0;
        sen_geo.auxil.e_sensor_layer[:,end] .= view(sun_geo.auxil.e_difꜛ,:,n_layer+1);
        sen_geo.auxil.e_sensor .= view(sun_geo.auxil.e_difꜛ,:,n_layer+1) ./ FT(π);
        sen_geo.auxil.reflectance .= sbulk.auxil.ρ_sw;

        return nothing
    end;

    # Run the canopy optical properties simulations only if canopy reflectance feature is enabled

    # compute the spectra at the observer direction
    for irt in 1:n_layer
        s_d_i = view(sun_geo.auxil.e_dirꜜ,:,irt);           # direct radiation at upper boundary
        e_d_i = view(sun_geo.auxil.e_difꜜ,:,irt);           # downward diffuse radiation at upper boundary
        e_u_j = view(sun_geo.auxil.e_difꜛ,:,irt+1);         # upward diffuse radiation at upper boundary
        sen_i = view(sen_geo.auxil.e_sensor_layer,:,irt);   # radiation towards the viewing direction per layer (including soil)

        #=
        ρ_do_layer = view(sen_geo.auxil.ρ_do_layer,:,irt);  # scattering coefficient from diffuse->observer
        τ_do_layer = view(sen_geo.auxil.τ_do_layer,:,irt);  # transmission coefficient from diffuse->observer
        ρ_so_layer = view(sen_geo.auxil.ρ_so_layer,:,irt);  # scattering coefficient from solar->observer

        Σlai = sum(view(can_str.trait.δlai,1:irt-1));
        Σsai = sum(view(can_str.trait.δsai,1:irt-1));
        kt_oo_x = sen_geo.auxil.ko_leaf * Σlai + sen_geo.auxil.ko_stem * Σsai;
        do_escape = exp(-kt_oo_x);
        do_escape = sen_geo.auxil.p_sensor[irt];
        sen_i .= do_escape .* (s_d_i .* ρ_so_layer .+ e_d_i .* ρ_do_layer) .+ sen_geo.auxil.p_sun_sensor[irt] ./ sun_geo.auxil.p_sunlit[irt] .* e_u_j .* τ_do_layer;
        =#



        # #=
        dob_l = view(sen_geo.auxil.dob_leaf,:,irt);       # scattering coefficient backward for diffuse->observer
        dof_l = view(sen_geo.auxil.dof_leaf,:,irt);       # scattering coefficient forward for diffuse->observer
        so_l  = view(sen_geo.auxil.so_leaf ,:,irt);       # bidirectional from solar to observer
        dob_s = view(sen_geo.auxil.dob_stem,:,irt);       # scattering coefficient backward for diffuse->observer
        dof_s = view(sen_geo.auxil.dof_stem,:,irt);       # scattering coefficient forward for diffuse->observer
        so_s  = view(sen_geo.auxil.so_stem ,:,irt);       # bidirectional from solar to observer

        # note here that ci is already accounted for in the p_sensor, so remove it from the equation here
        ilai = can_str.trait.δlai[irt];
        isai = can_str.trait.δsai[irt];
        sen_i .= sen_geo.auxil.p_sensor[irt] .* ilai .* sen_geo.auxil.ko_leaf .* (dob_l .* e_d_i .+ dof_l .* e_u_j) .+ sen_geo.auxil.p_sun_sensor[irt] .* ilai .* so_l .* rad_sw.e_dir .+
                 sen_geo.auxil.p_sensor[irt] .* isai .* sen_geo.auxil.ko_stem .* (dob_s .* e_d_i .+ dof_s .* e_u_j) .+ sen_geo.auxil.p_sun_sensor[irt] .* isai .* so_s .* rad_sw.e_dir;
        # =#
    end;



    sen_geo.auxil.e_sensor_layer[:,end] .= sen_geo.auxil.p_sensor_soil .* view(sun_geo.auxil.e_difꜛ,:,n_layer+1);

    # compute the spectra at the sensor
    for i in eachindex(sen_geo.auxil.e_sensor)
        sen_geo.auxil.e_sensor[i] = sum(view(sen_geo.auxil.e_sensor_layer,i,:)) / FT(π);
    end;

    # Note, this reflectance calculation is not correct because the sun-sensor geometry is not taken into account (reflectance is not isotropic)
    #     This is to compare with remote sensing data, which use the same calculation (* π)
    #     SCOPE does this a bit differently (controlling numerical issues), but still not correct
    # For real albedo, which is the surface reflectance, it needs to be the ratio between upward diffuse light and total downward light
    sen_geo.auxil.reflectance .= sen_geo.auxil.e_sensor .* FT(π) ./ (rad_sw.e_dir .+ rad_sw.e_dif);

    return nothing
end;

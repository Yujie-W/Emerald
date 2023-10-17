# This file contains functions to run shortwave radiation simulations

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-11: add function to compute shortwave radiation
#     2023-Oct-11: fix a typo when using fs_abs (was using fs_fo_abs)
#     2023-Oct-13: compute the canopy albedo using the integrated radiation of the entire hemisphere
#     2023-Oct-14: if LAI <= 0, run soil shortwave radiation only
#     2023-Oct-14: if SZA > 89, set all shortwave fluxes to 0
#
#######################################################################################################################################################################################################
"""

    shortwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update shortwave radiation related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function shortwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    (; CANOPY, LEAVES, METEO, SOIL_BULK) = spac;

    # if sza > 89, set all the radiation variables to 0
    if spac.CANOPY.sun_geometry.state.sza > 89
        CANOPY.sun_geometry.auxil.e_difꜜ .= 0;
        CANOPY.sun_geometry.auxil.e_difꜛ .= 0;
        CANOPY.sun_geometry.auxil.e_dirꜜ .= 0;
        CANOPY.sun_geometry.auxil.e_net_dif .= 0;
        CANOPY.sun_geometry.auxil.e_net_dir .= 0;
        CANOPY.sun_geometry.auxil.albedo .= NaN;
        CANOPY.sun_geometry.auxil.e_net_dif .= 0;
        CANOPY.sun_geometry.auxil.e_net_dir .= 0;
        CANOPY.sun_geometry.auxil.r_net_sw .= 0;
        SOIL_BULK.auxil.e_net_dir .= 0;
        SOIL_BULK.auxil.e_net_dif .= 0;
        SOIL_BULK.auxil.r_net_sw = 0;
        for leaf in LEAVES
            leaf.flux.auxil.apar_shaded = 0;
            leaf.flux.auxil.apar_sunlit .= 0;
            leaf.flux.auxil.ppar_shaded = 0;
            leaf.flux.auxil.ppar_sunlit .= 0;
        end;

        return nothing
    end;

    # if LAI <= 0, run soil albedo only
    (; DIM_LAYER, SPECTRA) = config;
    if spac.CANOPY.sun_geometry.state.sza <= 89 && spac.CANOPY.structure.state.lai <= 0
        # 1. update upward and downward direct and diffuse radiation profiles
        CANOPY.sun_geometry.auxil.e_dirꜜ .= METEO.rad_sw.e_dir;
        CANOPY.sun_geometry.auxil.e_difꜜ .= METEO.rad_sw.e_dif;
        CANOPY.sun_geometry.auxil.e_difꜛ .= 0;
        CANOPY.sun_geometry.auxil.e_difꜛ[:,end] .= (METEO.rad_sw.e_dir .+ METEO.rad_sw.e_dif) .* SOIL_BULK.auxil.ρ_sw;
        CANOPY.sun_geometry.auxil.albedo .= SOIL_BULK.auxil.ρ_sw;

        # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
        CANOPY.sun_geometry.auxil.e_net_dif .= 0;
        CANOPY.sun_geometry.auxil.e_net_dir .= 0;
        CANOPY.sun_geometry.auxil.r_net_sw .= 0;

        # 3. compute net absorption for leaves and soil
        SOIL_BULK.auxil.e_net_dir .= METEO.rad_sw.e_dir .* (1 .- SOIL_BULK.auxil.ρ_sw);
        SOIL_BULK.auxil.e_net_dif .= METEO.rad_sw.e_dif .* (1 .- SOIL_BULK.auxil.ρ_sw);
        SOIL_BULK.auxil.r_net_sw = (SOIL_BULK.auxil.e_net_dir' * SPECTRA.ΔΛ + SOIL_BULK.auxil.e_net_dif' * SPECTRA.ΔΛ) / 1000;

        # 4. compute leaf level PAR, APAR, and PPAR per ground area
        for leaf in LEAVES
            leaf.flux.auxil.apar_shaded = 0;
            leaf.flux.auxil.apar_sunlit .= 0;
            leaf.flux.auxil.ppar_shaded = 0;
            leaf.flux.auxil.ppar_sunlit .= 0;
        end;

        return nothing
    end;

    # run the shortwave radiation simulations only if LAI > 0 and when sza <= 89
    # 1. update upward and downward direct and diffuse radiation profiles
    CANOPY.sun_geometry.auxil.e_dirꜜ[:,1] .= METEO.rad_sw.e_dir;
    CANOPY.sun_geometry.auxil.e_difꜜ[:,1] .= METEO.rad_sw.e_dif;
    for i in 1:DIM_LAYER
        e_d_i = view(CANOPY.sun_geometry.auxil.e_difꜜ,:,i  );       # downward diffuse radiation at upper boundary
        e_d_j = view(CANOPY.sun_geometry.auxil.e_difꜜ,:,i+1);       # downward diffuse radiation at lower boundary
        e_s_i = view(CANOPY.sun_geometry.auxil.e_dirꜜ,:,i  );       # direct radiation at upper boundary
        e_s_j = view(CANOPY.sun_geometry.auxil.e_dirꜜ,:,i+1);       # direct radiation at lower boundary
        e_u_i = view(CANOPY.sun_geometry.auxil.e_difꜛ,:,i  );       # upward diffuse radiation at upper boundary

        r_dd_i = view(CANOPY.structure.auxil.ρ_dd         ,:,i);    # reflectance of the upper boundary (i)
        r_sd_i = view(CANOPY.sun_geometry.auxil.ρ_sd      ,:,i);    # reflectance of the upper boundary (i)
        t_dd_i = view(CANOPY.structure.auxil.τ_dd         ,:,i);    # transmittance of the layer (i)
        t_sd_i = view(CANOPY.sun_geometry.auxil.τ_sd      ,:,i);    # transmittance of the layer (i)
        t_ss_i = view(CANOPY.sun_geometry.auxil.τ_ss_layer,  i);    # transmittance for directional->directional

        e_s_j .= t_ss_i .* e_s_i;
        e_d_j .= t_sd_i .* e_s_i .+ t_dd_i .* e_d_i;
        e_u_i .= r_sd_i .* e_s_i .+ r_dd_i .* e_d_i;
    end;
    CANOPY.sun_geometry.auxil.e_difꜛ[:,end] .= view(CANOPY.sun_geometry.auxil.e_dirꜜ,:,DIM_LAYER+1) .* view(CANOPY.sun_geometry.auxil.ρ_sd,:,DIM_LAYER+1) .+
                                               view(CANOPY.sun_geometry.auxil.e_difꜜ,:,DIM_LAYER+1) .* view(CANOPY.structure.auxil.ρ_dd,:,DIM_LAYER+1);
    CANOPY.sun_geometry.auxil.albedo .= view(CANOPY.sun_geometry.auxil.e_difꜛ,:,1) ./ (METEO.rad_sw.e_dir .+ METEO.rad_sw.e_dif);

    # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
    for i in 1:DIM_LAYER
        e_d_i = view(CANOPY.sun_geometry.auxil.e_difꜜ,:,i  );       # downward diffuse radiation at upper boundary
        e_s_i = view(CANOPY.sun_geometry.auxil.e_dirꜜ,:,i  );       # direct radiation at upper boundary
        e_u_j = view(CANOPY.sun_geometry.auxil.e_difꜛ,:,i+1);       # upward diffuse radiation at upper boundary

        a_s_i = view(CANOPY.sun_geometry.auxil.e_net_dir,:,i);      # net absorbed direct radiation
        a_d_i = view(CANOPY.sun_geometry.auxil.e_net_dif,:,i);      # net absorbed diffuse radiation

        r_dd = view(CANOPY.structure.auxil.ρ_dd_layer   ,:,i);      # reflectance of the upper boundary (i)
        r_sd = view(CANOPY.sun_geometry.auxil.ρ_sd_layer,:,i);      # reflectance of the upper boundary (i)
        t_dd = view(CANOPY.structure.auxil.τ_dd_layer   ,:,i);      # transmittance of the layer (i)
        t_sd = view(CANOPY.sun_geometry.auxil.τ_sd_layer,:,i);      # transmittance of the layer (i)
        t_ss = view(CANOPY.sun_geometry.auxil.τ_ss_layer,  i);      # transmittance for directional->directional

        a_s_i .= e_s_i .* (1 .- t_ss .- t_sd .- r_sd);
        a_d_i .= (e_d_i .+ e_u_j) .* (1 .- t_dd .- r_dd);
    end;

    # 3. compute net absorption for leaves and soil
    for i in 1:DIM_LAYER
        Σ_shaded = view(CANOPY.sun_geometry.auxil.e_net_dif,:,i)' * SPECTRA.ΔΛ / 1000;
        Σ_sunlit = view(CANOPY.sun_geometry.auxil.e_net_dir,:,i)' * SPECTRA.ΔΛ / 1000;
        CANOPY.sun_geometry.auxil.r_net_sw[i] = (Σ_shaded + Σ_sunlit) / CANOPY.structure.state.δlai[i];
    end;
    SOIL_BULK.auxil.e_net_dir .= view(CANOPY.sun_geometry.auxil.e_dirꜜ,:,DIM_LAYER+1) .* (1 .- SOIL_BULK.auxil.ρ_sw);
    SOIL_BULK.auxil.e_net_dif .= view(CANOPY.sun_geometry.auxil.e_difꜜ,:,DIM_LAYER+1) .* (1 .- SOIL_BULK.auxil.ρ_sw);
    SOIL_BULK.auxil.r_net_sw = (SOIL_BULK.auxil.e_net_dir' * SPECTRA.ΔΛ + SOIL_BULK.auxil.e_net_dif' * SPECTRA.ΔΛ) / 1000;

    # 4. compute leaf level PAR, APAR, and PPAR per ground area
    normi = 1 / mean(CANOPY.sun_geometry.auxil.fs_abs_mean);
    for i in 1:DIM_LAYER
        j = DIM_LAYER + 1 - i;
        α_apar = view(LEAVES[j].bio.auxil.f_ppar, SPECTRA.IΛ_PAR);

        # convert energy to quantum unit for PAR, APAR and PPAR per leaf area
        CANOPY.sun_geometry.auxil._apar_shaded .= photon.(SPECTRA.Λ_PAR, view(CANOPY.sun_geometry.auxil.e_net_dif,SPECTRA.IΛ_PAR,i)) .* 1000 ./ CANOPY.structure.state.δlai[i];
        CANOPY.sun_geometry.auxil._apar_sunlit .= photon.(SPECTRA.Λ_PAR, view(CANOPY.sun_geometry.auxil.e_net_dir,SPECTRA.IΛ_PAR,i)) .* 1000 ./ CANOPY.structure.state.δlai[i] ./
                                                                                                                                                CANOPY.sun_geometry.auxil.p_sunlit[i];
        CANOPY.sun_geometry.auxil._ppar_shaded .= CANOPY.sun_geometry.auxil._apar_shaded .* α_apar;
        CANOPY.sun_geometry.auxil._ppar_sunlit .= CANOPY.sun_geometry.auxil._apar_sunlit .* α_apar;

        # APAR for leaves
        Σ_apar_dif = CANOPY.sun_geometry.auxil._apar_shaded' * SPECTRA.ΔΛ_PAR;
        Σ_apar_dir = CANOPY.sun_geometry.auxil._apar_sunlit' * SPECTRA.ΔΛ_PAR * normi;
        LEAVES[j].flux.auxil.apar_shaded = Σ_apar_dif;
        LEAVES[j].flux.auxil.apar_sunlit .= CANOPY.sun_geometry.auxil.fs_abs .* Σ_apar_dir .+ Σ_apar_dif;

        # PPAR for leaves
        Σ_ppar_dif = CANOPY.sun_geometry.auxil._ppar_shaded' * SPECTRA.ΔΛ_PAR;
        Σ_ppar_dir = CANOPY.sun_geometry.auxil._ppar_sunlit' * SPECTRA.ΔΛ_PAR * normi;
        LEAVES[j].flux.auxil.ppar_shaded = Σ_ppar_dif;
        LEAVES[j].flux.auxil.ppar_sunlit .= CANOPY.sun_geometry.auxil.fs_abs .* Σ_ppar_dir .+ Σ_ppar_dif;
    end;

    return nothing
end;
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
#     2023-Oct-18: account for SAI in the shortwave radiation calculation
#     2023-Oct-18: partition the energy between leaf and stem
#     2024-Jan-23: set PPAR to be the minimum of 2PPAR_700 and PPAR_750
#     2024-Jun-06: add step to compute e_difꜛ_layer (contribution from the layer only)
#     2024-Jul-27: use bined PPAR to speed up
#     2024-Jul-30: do not bin PPAR if DIM_PPAR_BINS is nothing
#     2024-Aug-05: add a special case when DIM_PPAR_BINS is 0 (one leaf model, no sunlit and shaded fraction)
#
#######################################################################################################################################################################################################
"""

    shortwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update shortwave radiation related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function shortwave_radiation! end;

shortwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = shortwave_radiation!(config, spac, spac.plant.leaves[1]);

shortwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, ::CanopyLayer{FT}) where {FT} = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    sbulk = spac.soil_bulk;
    sun_geo = spac.canopy.sun_geometry;
    n_layer = length(leaves);

    # if sza > 89, set all the radiation variables to 0
    if sun_geo.state.sza > 89
        sun_geo.auxil.e_difꜜ .= 0;
        sun_geo.auxil.e_difꜛ .= 0;
        sun_geo.auxil.e_dirꜜ .= 0;
        sun_geo.auxil.e_net_dif .= 0;
        sun_geo.auxil.e_net_dir .= 0;
        sun_geo.auxil.albedo .= NaN;
        sun_geo.auxil.e_net_dif .= 0;
        sun_geo.auxil.e_net_dir .= 0;
        sun_geo.auxil.r_net_sw_leaf .= 0;
        sun_geo.auxil.r_net_sw_stem .= 0;
        sbulk.auxil.e_net_dir .= 0;
        sbulk.auxil.e_net_dif .= 0;
        sbulk.auxil.r_net_sw = 0;
        for leaf in leaves
            leaf.flux.auxil.ppar .= 0;
        end;

        return nothing
    end;

    # if LAI <= 0, run soil albedo only
    (; DIM_AZI, DIM_INCL, DIM_PPAR_BINS, SPECTRA) = config;
    rad_sw = spac.meteo.rad_sw;
    if can_str.trait.lai <= 0 && can_str.trait.sai <= 0
        # 1. update upward and downward direct and diffuse radiation profiles
        sun_geo.auxil.e_dirꜜ .= rad_sw.e_dir;
        sun_geo.auxil.e_difꜜ .= rad_sw.e_dif;
        sun_geo.auxil.e_difꜛ .= 0;
        sun_geo.auxil.e_difꜛ[:,end] .= (rad_sw.e_dir .+ rad_sw.e_dif) .* sbulk.auxil.ρ_sw;
        sun_geo.auxil.albedo .= sbulk.auxil.ρ_sw;

        # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
        sun_geo.auxil.e_net_dif .= 0;
        sun_geo.auxil.e_net_dir .= 0;
        sun_geo.auxil.r_net_sw_leaf .= 0;
        sun_geo.auxil.r_net_sw_stem .= 0;

        # 3. compute net absorption for leaves and soil
        sbulk.auxil.e_net_dir .= rad_sw.e_dir .* (1 .- sbulk.auxil.ρ_sw);
        sbulk.auxil.e_net_dif .= rad_sw.e_dif .* (1 .- sbulk.auxil.ρ_sw);
        sbulk.auxil.r_net_sw = (sbulk.auxil.e_net_dir' * SPECTRA.ΔΛ + sbulk.auxil.e_net_dif' * SPECTRA.ΔΛ) / 1000;

        # 4. compute leaf level PAR, APAR, and PPAR per ground area
        for leaf in leaves
            leaf.flux.auxil.ppar .= 0;
        end;

        return nothing
    end;

    # run the shortwave radiation simulations only if LAI > 0 and when sza <= 89
    # 1. update upward and downward direct and diffuse radiation profiles
    sun_geo.auxil.e_dirꜜ[:,1] .= rad_sw.e_dir;
    sun_geo.auxil.e_difꜜ[:,1] .= rad_sw.e_dif;
    for i in 1:n_layer
        e_d_i = view(sun_geo.auxil.e_difꜜ      ,:,i  );    # downward diffuse radiation at upper boundary
        e_d_j = view(sun_geo.auxil.e_difꜜ      ,:,i+1);    # downward diffuse radiation at lower boundary
        e_s_i = view(sun_geo.auxil.e_dirꜜ      ,:,i  );    # direct radiation at upper boundary
        e_s_j = view(sun_geo.auxil.e_dirꜜ      ,:,i+1);    # direct radiation at lower boundary
        e_u_i = view(sun_geo.auxil.e_difꜛ      ,:,i  );    # upward diffuse radiation at upper boundary
        e_a_i = view(sun_geo.auxil.e_difꜛ_layer,:,i  );    # upward diffuse radiation at upper boundary (contribution from the layer only)

        r_dd_i = view(can_str.auxil.ρ_dd      ,:,i  );  # reflectance of the upper boundary (i)
        r_sd_i = view(sun_geo.auxil.ρ_sd      ,:,i  );  # reflectance of the upper boundary (i)
        t_dd_i = view(can_str.auxil.τ_dd      ,:,i  );  # transmittance of the layer (i)
        t_sd_i = view(sun_geo.auxil.τ_sd      ,:,i  );  # transmittance of the layer (i)
        t_ss__ = view(sun_geo.auxil.τ_ss_layer,  i  );  # transmittance for directional->directional
        r_sd__ = view(sun_geo.auxil.ρ_sd_layer,:,i  );  # reflectance for directional->diffuse
        r_dd__ = view(can_str.auxil.ρ_dd_layer,:,i  );  # reflectance for diffuse->diffuse

        e_s_j .= t_ss__ .* e_s_i;
        e_d_j .= t_sd_i .* e_s_i .+ t_dd_i .* e_d_i;
        e_u_i .= r_sd_i .* e_s_i .+ r_dd_i .* e_d_i;
        e_a_i .= r_sd__ .* e_s_i .+ r_dd__ .* e_d_i;
    end;
    sun_geo.auxil.e_difꜛ[:,end] .= view(sun_geo.auxil.e_dirꜜ,:,n_layer+1) .* view(sun_geo.auxil.ρ_sd,:,n_layer+1) .+
                                   view(sun_geo.auxil.e_difꜜ,:,n_layer+1) .* view(can_str.auxil.ρ_dd,:,n_layer+1);
    sun_geo.auxil.e_difꜛ_layer[:,end] .= view(sun_geo.auxil.e_difꜛ,:,size(sun_geo.auxil.e_difꜛ,2));
    sun_geo.auxil.albedo .= view(sun_geo.auxil.e_difꜛ,:,1) ./ (rad_sw.e_dir .+ rad_sw.e_dif);

    # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
    for i in 1:n_layer
        e_d_i = view(sun_geo.auxil.e_difꜜ,:,i  );       # downward diffuse radiation at upper boundary
        e_s_i = view(sun_geo.auxil.e_dirꜜ,:,i  );       # direct radiation at upper boundary
        e_u_j = view(sun_geo.auxil.e_difꜛ,:,i+1);       # upward diffuse radiation at upper boundary

        a_s_i = view(sun_geo.auxil.e_net_dir,:,i);      # net absorbed direct radiation
        a_d_i = view(sun_geo.auxil.e_net_dif,:,i);      # net absorbed diffuse radiation

        r_dd = view(can_str.auxil.ρ_dd_layer,:,i);      # reflectance of the upper boundary (i)
        r_sd = view(sun_geo.auxil.ρ_sd_layer,:,i);      # reflectance of the upper boundary (i)
        t_dd = view(can_str.auxil.τ_dd_layer,:,i);      # transmittance of the layer (i)
        t_sd = view(sun_geo.auxil.τ_sd_layer,:,i);      # transmittance of the layer (i)
        t_ss = view(sun_geo.auxil.τ_ss_layer,  i);      # transmittance for directional->directional

        a_s_i .= e_s_i .* (1 .- t_ss .- t_sd .- r_sd);
        a_d_i .= (e_d_i .+ e_u_j) .* (1 .- t_dd .- r_dd);
    end;

    # 3. compute net absorption for leaves and soil
    # 4. compute leaf level PAR, APAR, and PPAR per ground area
    normi = 1 / mean(sun_geo.s_aux.fs_abs_mean);
    a_leaf = spac.cache.cache_wl_1;
    a_stem = spac.cache.cache_wl_2;
    f_leaf = spac.cache.cache_wl_3;
    f_stem = spac.cache.cache_wl_4;
    temp_p = spac.cache.cache_wl_5;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        a_leaf .= leaf.bio.auxil.α_leaf .* can_str.trait.δlai[irt];
        a_stem .= (1 .- SPECTRA.ρ_STEM) .* can_str.trait.δsai[irt];
        f_leaf .= a_leaf ./ (a_leaf .+ a_stem);
        f_stem .= 1 .- f_leaf;

        temp_p .= view(sun_geo.auxil.e_net_dif,:,irt) .* f_leaf;
        Σ_shaded_leaf = temp_p' * SPECTRA.ΔΛ / 1000;
        temp_p .= view(sun_geo.auxil.e_net_dir,:,irt) .* f_leaf;
        Σ_sunlit_leaf = temp_p' * SPECTRA.ΔΛ / 1000;
        temp_p .= view(sun_geo.auxil.e_net_dif,:,irt) .* f_stem;
        Σ_shaded_stem = temp_p' * SPECTRA.ΔΛ / 1000;
        temp_p .= view(sun_geo.auxil.e_net_dir,:,irt) .* f_stem;
        Σ_sunlit_stem = temp_p' * SPECTRA.ΔΛ / 1000;

        # partition the net radiation to leaves and stem
        sun_geo.auxil.r_net_sw_leaf[irt] = Σ_shaded_leaf + Σ_sunlit_leaf;
        sun_geo.auxil.r_net_sw_stem[irt] = Σ_shaded_stem + Σ_sunlit_stem;

        # compute leaf level PAR, APAR, and PPAR per ground area
        if can_str.trait.δlai[irt] > 0
            α_apar = view(leaf.bio.auxil.f_ppar, SPECTRA.IΛ_PAR);
            p_leaf = view(f_leaf, SPECTRA.IΛ_PAR);
            # convert energy to quantum unit for PAR, APAR and PPAR per leaf area
            sun_geo.auxil._apar_shaded .= photon.(SPECTRA.Λ_PAR, view(sun_geo.auxil.e_net_dif,SPECTRA.IΛ_PAR,irt)) .* p_leaf .* 1000 ./ can_str.trait.δlai[irt];
            sun_geo.auxil._apar_sunlit .= photon.(SPECTRA.Λ_PAR, view(sun_geo.auxil.e_net_dir,SPECTRA.IΛ_PAR,irt)) .* p_leaf .* 1000 ./ can_str.trait.δlai[irt] ./ sun_geo.s_aux.p_sunlit[irt];
            sun_geo.auxil._ppar_shaded .= sun_geo.auxil._apar_shaded .* α_apar;
            sun_geo.auxil._ppar_sunlit .= sun_geo.auxil._apar_sunlit .* α_apar;

            # APAR for leaves
            Σ_apar_dif = sun_geo.auxil._apar_shaded' * SPECTRA.ΔΛ_PAR;
            Σ_apar_dir = sun_geo.auxil._apar_sunlit' * SPECTRA.ΔΛ_PAR * normi;
            sun_geo.auxil.apar_shaded[irt] = Σ_apar_dif;
            sun_geo.auxil.apar_sunlit[:,:,irt] .= sun_geo.s_aux.fs_abs .* Σ_apar_dir .+ Σ_apar_dif;

            # PPAR for leaves (set PPAR to be the minimum of 2PPAR_700 and PPAR_750)
            Σ_ppar_dif_700 = view(sun_geo.auxil._ppar_shaded, SPECTRA.IΛ_PAR_700)' * SPECTRA.ΔΛ_PAR_700;
            Σ_ppar_dir_700 = view(sun_geo.auxil._ppar_sunlit, SPECTRA.IΛ_PAR_700)' * SPECTRA.ΔΛ_PAR_700 * normi;
            Σ_ppar_dif_750 = sun_geo.auxil._ppar_shaded' * SPECTRA.ΔΛ_PAR;
            Σ_ppar_dir_750 = sun_geo.auxil._ppar_sunlit' * SPECTRA.ΔΛ_PAR * normi;
            Σ_ppar_dif = min(2Σ_ppar_dif_700, Σ_ppar_dif_750);
            Σ_ppar_dir = min(2Σ_ppar_dir_700, Σ_ppar_dir_750);
            sun_geo.auxil.ppar_shaded[irt] = Σ_ppar_dif;
            sun_geo.auxil.ppar_sunlit[:,:,irt] .= sun_geo.s_aux.fs_abs .* Σ_ppar_dir .+ Σ_ppar_dif;

            # bin the PPAR values based on their PPAR if DIM_PPAR_BINS is not nothing
            if isnothing(DIM_PPAR_BINS)
                for j in 1:DIM_AZI
                    leaf.flux.auxil.ppar[(j-1)*DIM_INCL+1:j*DIM_INCL] .= view(sun_geo.s_aux.fs_abs,:,j) .* Σ_ppar_dir .+ Σ_ppar_dif;
                end;
                sun_geo.auxil.ppar_fraction[1:end-1,irt] .= 1 ./ ( DIM_INCL * DIM_AZI) .* sun_geo.s_aux.p_sunlit[irt];
            elseif DIM_PPAR_BINS == 0
                # mean_ppar_sunlit = mean(view(sun_geo.auxil.ppar_sunlit,:,:,irt));
                # leaf.flux.auxil.ppar[1] = mean_ppar_sunlit * sun_geo.s_aux.p_sunlit[irt] + Σ_ppar_dif * (1 - sun_geo.s_aux.p_sunlit[irt]);
                leaf.flux.auxil.ppar[1] = Σ_ppar_dir / normi * sun_geo.s_aux.p_sunlit[irt] + Σ_ppar_dif;
                sun_geo.auxil.ppar_index .= 1;
            else
                min_ppar = minimum(view(sun_geo.auxil.ppar_sunlit, :, :, irt));
                max_ppar = maximum(view(sun_geo.auxil.ppar_sunlit, :, :, irt));
                δ_ppar = (max_ppar - min_ppar) / DIM_PPAR_BINS;
                sun_geo.auxil._ppar_sum .= 0;
                sun_geo.auxil._ppar_count .= 0;
                if δ_ppar == 0
                    sun_geo.auxil._ppar_sum[1] += max_ppar * DIM_INCL * DIM_AZI;
                    sun_geo.auxil._ppar_count[1] += DIM_INCL * DIM_AZI;
                else
                    for i in 1:DIM_INCL, j in 1:DIM_AZI
                        ind = min(DIM_PPAR_BINS, Int((sun_geo.auxil.ppar_sunlit[i,j,irt] - min_ppar) ÷ δ_ppar + 1));
                        sun_geo.auxil._ppar_sum[ind] += sun_geo.auxil.ppar_sunlit[i,j,irt];
                        sun_geo.auxil._ppar_count[ind] += 1;
                        sun_geo.auxil.ppar_index[i,j,irt] = ind;
                    end;
                end;
                for i in 1:DIM_PPAR_BINS
                    if sun_geo.auxil._ppar_count[i] > 0
                        leaf.flux.auxil.ppar[i] = sun_geo.auxil._ppar_sum[i] / sun_geo.auxil._ppar_count[i];
                    else
                        leaf.flux.auxil.ppar[i] = 0;
                    end;
                end;
                sun_geo.auxil.ppar_fraction[1:DIM_PPAR_BINS,irt] .= sun_geo.auxil._ppar_count ./ ( DIM_INCL * DIM_AZI) .* sun_geo.s_aux.p_sunlit[irt];
            end;
            leaf.flux.auxil.ppar[end] = Σ_ppar_dif;
            sun_geo.auxil.ppar_fraction[end,irt] = 1 - sun_geo.s_aux.p_sunlit[irt];
        else
            leaf.flux.auxil.ppar .= 0;
        end;
    end;
    sbulk.auxil.e_net_dir .= view(sun_geo.auxil.e_dirꜜ,:,n_layer+1) .* (1 .- sbulk.auxil.ρ_sw);
    sbulk.auxil.e_net_dif .= view(sun_geo.auxil.e_difꜜ,:,n_layer+1) .* (1 .- sbulk.auxil.ρ_sw);
    sbulk.auxil.r_net_sw = (sbulk.auxil.e_net_dir' * SPECTRA.ΔΛ + sbulk.auxil.e_net_dif' * SPECTRA.ΔΛ) / 1000;

    return nothing
);

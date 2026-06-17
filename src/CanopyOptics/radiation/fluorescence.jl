# This file contains functions to compute the SIF emission of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-14: add function fluorescence_spectrum! (run per sensor geometry)
#     2023-Oct-14: if LAI < = 0 or SZA > 89, set all fluxes to 0
#     2023-Oct-18: SIF excitation is rescaled to leaf partitioning (accounting stem)
#     2024-Jun-07: add step to compute e_sifꜛ_layer_sum (contribution to upward SIF from the layer after relection from the lower layers)
#     2024-Jul-27: use bined PPAR to speed up
#     2024-Jul-30: do not bin PPAR if DIM_PPAR_BINS is nothing
#     2024-Sep-04: separate leaf and stem optical properties
# Bug fixes
#     2024-Mar-06: ci impact on fraction from viewer direction (otherwise will be accounted twice)
#     2024-Sep-07: do not use CI in the SIF emission calculation (introduced in 2024-Mar-06)
#     2025-Sep-12: add a special case when toral rad is zero (to avoid NaN issue)
#     2026-Mar-27: use e_difꜛ from jth layer for fluoresence excitation of the ith layer (used the e_difꜛ from ith layer, which is incorrect)
#
#######################################################################################################################################################################################################
"""

    fluorescence_spectrum!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Compute the fluorescence spectrum of the canopy at the sensor direction, given
- `config` SPAC configuration
- `spac` SPAC

"""
function fluorescence_spectrum!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    if !config.FEATURES.ENABLE_SIF
        return nothing
    end;

    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    sun_geo = spac.canopy.sun_geometry;
    n_layer = length(leaves);
    rad_sw = spac.meteo.rad_sw;
    (; SPECTRA) = config.CONSTANTS;
    (; DIM_AZI, DIM_INCL, DIM_PPAR_BINS) = config.DIMENSIONS;

    # if sza > 89, set all the radiation variables to 0
    total_sw_rad = (rad_sw.e_dir' * SPECTRA.ΔΛ + rad_sw.e_dif' * SPECTRA.ΔΛ) / 1000;
    if sun_geo.state.sza > 89 || can_str.trait.lai <= 0 || total_sw_rad <= 0
        sun_geo.auxil.e_sif_chl .= 0;
        sun_geo.auxil.e_sifꜜ_layer .= 0;
        sun_geo.auxil.e_sifꜛ_layer .= 0;
        sun_geo.auxil.e_sifꜜ_emit .= 0;
        sun_geo.auxil.e_sifꜛ_emit .= 0;
        sun_geo.auxil.e_sifꜜ .= 0;
        sun_geo.auxil.e_sifꜛ .= 0;
        sen_geo.auxil.sif_sunlit .= 0;
        sen_geo.auxil.sif_shaded .= 0;
        sen_geo.auxil.sif_scattered .= 0;
        sen_geo.auxil.sif_obs_sunlit .= 0;
        sen_geo.auxil.sif_obs_shaded .= 0;
        sen_geo.auxil.sif_obs_scattered .= 0;
        sen_geo.auxil.sif_obs_soil .= 0;
        sen_geo.auxil.sif_obs .= 0;

        return nothing
    end;

    # run the fluorescence simulations only if fluorescence feature is enabled
    # broadcast the phi_f from bined array to 3D array
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sen_geo.auxil.ϕ_f_shaded[irt] = leaf.photosystem.auxil.ϕ_f[end];
        if isnothing(DIM_PPAR_BINS)
            for i in 1:DIM_AZI
                sen_geo.auxil.ϕ_f_sunlit[irt][:,i] .= view(leaf.photosystem.auxil.ϕ_f,(i-1)*DIM_INCL+1:i*DIM_INCL);
            end;
        else
            for i in 1:DIM_INCL, j in 1:DIM_AZI
                sen_geo.auxil.ϕ_f_sunlit[irt][i,j] = leaf.photosystem.auxil.ϕ_f[ sun_geo.auxil.ppar_index[i,j,irt] ];
            end;
        end;
    end;

    # function to weight matrices by inclination angles
    @inline local_lidf_weight(mat_0, mat_1) = (
        sun_geo.auxil._mat_incl_azi .= mat_0 .* mat_1;
        mul!(sun_geo.auxil._vec_azi, sun_geo.auxil._mat_incl_azi', can_str.auxil.p_incl_leaf);

        return mean(sun_geo.auxil._vec_azi)
    );
    @inline local_lidf_weight(mat_1) = (
        mul!(sun_geo.auxil._vec_azi, mat_1', can_str.auxil.p_incl_leaf);

        return mean(sun_geo.auxil._vec_azi)
    );
    _COS²_Θ_INCL_AZI = spac.cache.cache_incl_azi_1;
    _COS²_Θ_INCL_AZI .= (cosd.(config.DIMENSIONS.Θ_INCL) .^ 2);

    #
    #
    # TODO: better use mat_b and mat_f for a non-reabsorbing scenario
    #
    #
    # 0. compute chloroplast SIF emissions for different layers
    a_leaf = spac.cache.cache_sife_1;
    a_stem = spac.cache.cache_sife_2;
    f_leaf = spac.cache.cache_sife_3;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        a_leaf .= view(leaf.bio.auxil.α_leaf,SPECTRA.IΛ_SIFE) .* can_str.trait.δlai[irt];
        a_stem .= (1 .- view(SPECTRA.ρ_STEM,SPECTRA.IΛ_SIFE)) .* can_str.trait.δsai[irt];
        f_leaf .= a_leaf ./ (a_leaf .+ a_stem);

        # compute the energy used for SIF excitation
        sun_geo.auxil._e_dirꜜ_sife .= view(sun_geo.auxil.e_dirꜜ,SPECTRA.IΛ_SIFE,irt  ) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜜ_sife .= view(sun_geo.auxil.e_difꜜ,SPECTRA.IΛ_SIFE,irt  ) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜛ_sife .= view(sun_geo.auxil.e_difꜛ,SPECTRA.IΛ_SIFE,irt+1) .* f_leaf .* SPECTRA.ΔΛ_SIFE;

        # convert the excitation radiation to photons if ϕ_photon is true
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_dirꜜ_sife);
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜜ_sife);
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜛ_sife);

        # convert the excitation radiation to fluorescence components
        mul!(sun_geo.auxil._e_dirꜜ_sifꜛ, leaf.bio.auxil.matꜛ_chl, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_dirꜜ_sifꜜ, leaf.bio.auxil.matꜜ_chl, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sifꜛ, leaf.bio.auxil.matꜛ_chl, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sifꜜ, leaf.bio.auxil.matꜜ_chl, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sifꜛ, leaf.bio.auxil.matꜛ_chl, sun_geo.auxil._e_difꜛ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sifꜜ, leaf.bio.auxil.matꜜ_chl, sun_geo.auxil._e_difꜛ_sife);

        # convert the SIF back to energy unit if ϕ_photon is true
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sifꜜ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sifꜜ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sifꜜ);

        #
        #
        # TODO: refactor this part when fully understand what is happening here
        #
        #
        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        ϕ_sunlit = sen_geo.auxil.ϕ_f_sunlit[irt];
        ϕ_shaded = sen_geo.auxil.ϕ_f_shaded[irt];

        # compute the weights
        sl_1_ = local_lidf_weight(ϕ_sunlit, 1);                             # SCOPE: etau_lidf
        sh_1_ = local_lidf_weight(ϕ_shaded, 1);                             # SCOPE: etah_lidf
        sl_θ² = local_lidf_weight(ϕ_sunlit, _COS²_Θ_INCL_AZI);              # SCOPE: bsxfun(@times,etau_lidf,ctl2)
        sh_θ² = local_lidf_weight(ϕ_shaded, _COS²_Θ_INCL_AZI);              # SCOPE: bsxfun(@times,etah_lidf,ctl2)
        sl_S_ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_abs);          # SCOPE: bsxfun(@times,etau_lidf,absfs)
        sl_sθ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_cos_incl);     # SCOPE: bsxfun(@times,etau_lidf,fsctl)

        # upward and downward SIF from direct and diffuse radiation per leaf area
        sun_geo.auxil._sif_sunlitꜛ_dif .= sun_geo.auxil._e_difꜜ_sifꜛ .* sl_1_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sl_θ² .+       # SCOPE: sigbEmin_u
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sl_1_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sl_θ²;         # SCOPE: sigfEplu_u
        sun_geo.auxil._sif_sunlitꜜ_dif .= sun_geo.auxil._e_difꜜ_sifꜛ .* sl_1_ .- sun_geo.auxil._e_difꜜ_sifꜜ .* sl_θ² .+       # SCOPE: sigfEmin_u
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sl_1_ .+ sun_geo.auxil._e_difꜛ_sifꜜ .* sl_θ²;         # SCOPE: sigbEplu_u
        sun_geo.auxil._sif_sunlitꜛ_dir .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_S_ .+ sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_sθ;         # SCOPE: sbEs
        sun_geo.auxil._sif_sunlitꜜ_dir .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_S_ .- sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_sθ;         # SCOPE: sfEs
        sun_geo.auxil._sif_shadedꜛ     .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_1_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sh_θ² .+       # SCOPE: sigbEmin_h
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sh_1_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sh_θ²;         # SCOPE: sigfEplu_h
        sun_geo.auxil._sif_shadedꜜ     .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_1_ .- sun_geo.auxil._e_difꜜ_sifꜜ .* sh_θ² .+       # SCOPE: sigfEmin_h
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sh_1_ .+ sun_geo.auxil._e_difꜛ_sifꜜ .* sh_θ²;         # SCOPE: sigbEplu_h

        # total emitted SIF for upward and downward direction (ci is already accounted for in p_sunlit, p_sun_sensor, and shortwave radiation, and thus there is no need to use CI here)
        # add ci_diffuse back to account for the scattering within the canopy layer
        # TODO: better SIF scattering algorithm
        # ciilai = (1 - exp(-can_str.trait.δlai[irt])) * can_str.auxil.ci_diffuse;
        # ciilai = can_str.trait.δlai[irt] * can_str.auxil.ci_diffuse;
        ilai_direct = (1 - sun_geo.auxil.τ_ss_layer[irt]) / sun_geo.auxil.ks_leaf;
        ilai_diffuse = 1 - can_str.auxil.τ_dd_isotropic[irt];
        sun_geo.auxil.e_sifꜜ_layer[:,irt] .= sun_geo.auxil._sif_sunlitꜜ_dir .* ilai_direct .+
                                             sun_geo.auxil._sif_sunlitꜜ_dif .* sun_geo.auxil.p_sunlit[irt] .* ilai_diffuse .+
                                             sun_geo.auxil._sif_shadedꜜ     .* (1 - sun_geo.auxil.p_sunlit[irt]) .* ilai_diffuse;
        sun_geo.auxil.e_sifꜛ_layer[:,irt] .= sun_geo.auxil._sif_sunlitꜛ_dir .* ilai_direct .+
                                             sun_geo.auxil._sif_sunlitꜛ_dif .* sun_geo.auxil.p_sunlit[irt] .* ilai_diffuse .+
                                             sun_geo.auxil._sif_shadedꜛ     .* (1 - sun_geo.auxil.p_sunlit[irt]) .* ilai_diffuse;
    end;
    sun_geo.auxil.e_sif_chl .= sun_geo.auxil.e_sifꜜ_layer .+ sun_geo.auxil.e_sifꜛ_layer;

    # 1. compute SIF emissions for different layers
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        a_leaf .= view(leaf.bio.auxil.α_leaf,SPECTRA.IΛ_SIFE) .* can_str.trait.δlai[irt];
        a_stem .= (1 .- view(SPECTRA.ρ_STEM,SPECTRA.IΛ_SIFE)) .* can_str.trait.δsai[irt];
        f_leaf .= a_leaf ./ (a_leaf .+ a_stem);

        # compute the energy used for SIF excitation
        sun_geo.auxil._e_dirꜜ_sife .= view(sun_geo.auxil.e_dirꜜ,SPECTRA.IΛ_SIFE,irt  ) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜜ_sife .= view(sun_geo.auxil.e_difꜜ,SPECTRA.IΛ_SIFE,irt  ) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜛ_sife .= view(sun_geo.auxil.e_difꜛ,SPECTRA.IΛ_SIFE,irt+1) .* f_leaf .* SPECTRA.ΔΛ_SIFE;

        # convert the excitation radiation to photons if ϕ_photon is true
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_dirꜜ_sife);
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜜ_sife);
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜛ_sife);

        # convert the excitation radiation to fluorescence components
        mul!(sun_geo.auxil._e_dirꜜ_sifꜛ, leaf.bio.auxil.matꜛ, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_dirꜜ_sifꜜ, leaf.bio.auxil.matꜜ, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sifꜛ, leaf.bio.auxil.matꜛ, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sifꜜ, leaf.bio.auxil.matꜜ, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sifꜛ, leaf.bio.auxil.matꜛ, sun_geo.auxil._e_difꜛ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sifꜜ, leaf.bio.auxil.matꜜ, sun_geo.auxil._e_difꜛ_sife);

        # convert the SIF back to energy unit if ϕ_photon is true
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sifꜜ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sifꜜ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sifꜜ);

        #
        #
        # TODO: refactor this part when fully understand what is happening here
        #
        #
        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        ϕ_sunlit = sen_geo.auxil.ϕ_f_sunlit[irt];
        ϕ_shaded = sen_geo.auxil.ϕ_f_shaded[irt];

        # compute the weights
        sl_1_ = local_lidf_weight(ϕ_sunlit, 1);                             # SCOPE: etau_lidf
        sh_1_ = local_lidf_weight(ϕ_shaded, 1);                             # SCOPE: etah_lidf
        sl_θ² = local_lidf_weight(ϕ_sunlit, _COS²_Θ_INCL_AZI);              # SCOPE: bsxfun(@times,etau_lidf,ctl2)
        sh_θ² = local_lidf_weight(ϕ_shaded, _COS²_Θ_INCL_AZI);              # SCOPE: bsxfun(@times,etah_lidf,ctl2)
        sl_S_ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_abs);          # SCOPE: bsxfun(@times,etau_lidf,absfs)
        sl_sθ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_cos_incl);     # SCOPE: bsxfun(@times,etau_lidf,fsctl)

        sh_O_ = local_lidf_weight(ϕ_shaded, sen_geo.auxil.fo_abs);          # SCOPE: bsxfun(@times,etah_lidf,absfo)
        sl_O_ = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_abs);          # SCOPE: bsxfun(@times,etau_lidf,absfo)
        sh_oθ = local_lidf_weight(ϕ_shaded, sen_geo.auxil.fo_cos_incl);     # SCOPE: bsxfun(@times,etah_lidf,foctl)
        sl_oθ = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_cos_incl);     # SCOPE: bsxfun(@times,etau_lidf,foctl)
        sl_SO = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_fs_abs);       # SCOPE: bsxfun(@times,etau_lidf,absfsfo)
        sl_so = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_fs);           # SCOPE: bsxfun(@times,etau_lidf,fsfo)

        # upward and downward SIF from direct and diffuse radiation per leaf area
        sun_geo.auxil._sif_sunlitꜛ_dif .= sun_geo.auxil._e_difꜜ_sifꜛ .* sl_1_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sl_θ² .+       # SCOPE: sigbEmin_u
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sl_1_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sl_θ²;         # SCOPE: sigfEplu_u
        sun_geo.auxil._sif_sunlitꜜ_dif .= sun_geo.auxil._e_difꜜ_sifꜛ .* sl_1_ .- sun_geo.auxil._e_difꜜ_sifꜜ .* sl_θ² .+       # SCOPE: sigfEmin_u
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sl_1_ .+ sun_geo.auxil._e_difꜛ_sifꜜ .* sl_θ²;         # SCOPE: sigbEplu_u
        sun_geo.auxil._sif_sunlitꜛ_dir .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_S_ .+ sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_sθ;         # SCOPE: sbEs
        sun_geo.auxil._sif_sunlitꜜ_dir .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_S_ .- sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_sθ;         # SCOPE: sfEs
        sun_geo.auxil._sif_shadedꜛ     .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_1_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sh_θ² .+       # SCOPE: sigbEmin_h
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sh_1_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sh_θ²;         # SCOPE: sigfEplu_h
        sun_geo.auxil._sif_shadedꜜ     .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_1_ .- sun_geo.auxil._e_difꜜ_sifꜜ .* sh_θ² .+       # SCOPE: sigfEmin_h
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sh_1_ .+ sun_geo.auxil._e_difꜛ_sifꜜ .* sh_θ²;         # SCOPE: sigbEplu_h

        # update the SIF cache for the observer direction (compute it here to save time)
        sen_geo.auxil.sif_sunlit[:,irt] .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_SO .+ sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_so .+      # SCOPE: wfEs
                                           sun_geo.auxil._e_difꜜ_sifꜛ .* sl_O_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sl_oθ .+      # SCOPE: vbEmin_u
                                           sun_geo.auxil._e_difꜛ_sifꜛ .* sl_O_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sl_oθ;        # SCOPE: vfEplu_u
        sen_geo.auxil.sif_shaded[:,irt] .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_O_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sh_oθ .+      # SCOPE: vbEmin_h
                                           sun_geo.auxil._e_difꜛ_sifꜛ .* sh_O_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sh_oθ;        # SCOPE: vfEplu_h

        # total emitted SIF for upward and downward direction (ci is already accounted for in p_sunlit, p_sun_sensor, and shortwave radiation, and thus there is no need to use CI here)
        # add ci_diffuse back to account for the scattering within the canopy layer
        # TODO: better SIF scattering algorithm
        # ciilai = (1 - exp(-can_str.trait.δlai[irt])) * can_str.auxil.ci_diffuse;
        # ciilai = can_str.trait.δlai[irt] * can_str.auxil.ci_diffuse;
        ilai_direct = (1 - sun_geo.auxil.τ_ss_layer[irt]) / sun_geo.auxil.ks_leaf;
        ilai_diffuse = 1 - can_str.auxil.τ_dd_isotropic[irt];
        sun_geo.auxil.e_sifꜜ_layer[:,irt] .= sun_geo.auxil._sif_sunlitꜜ_dir .* ilai_direct .+
                                             sun_geo.auxil._sif_sunlitꜜ_dif .* sun_geo.auxil.p_sunlit[irt] .* ilai_diffuse .+
                                             sun_geo.auxil._sif_shadedꜜ     .* (1 - sun_geo.auxil.p_sunlit[irt]) .* ilai_diffuse;
        sun_geo.auxil.e_sifꜛ_layer[:,irt] .= sun_geo.auxil._sif_sunlitꜛ_dir .* ilai_direct .+
                                             sun_geo.auxil._sif_sunlitꜛ_dif .* sun_geo.auxil.p_sunlit[irt] .* ilai_diffuse .+
                                             sun_geo.auxil._sif_shadedꜛ     .* (1 - sun_geo.auxil.p_sunlit[irt]) .* ilai_diffuse;
    end;

    # 2. account for the SIF emission from bottom to up
    sun_geo.auxil.e_sifꜛ_emit[:,end] .= 0;
    for i in n_layer:-1:1
        r__ = view(can_str.auxil.ρ_dd_layer,SPECTRA.IΛ_SIF,i  );    # reflectance of the layer without correction
        r_j = view(can_str.auxil.ρ_dd      ,SPECTRA.IΛ_SIF,i+1);    # reflectance of the lower boundary (i) for SIF
        t__ = view(can_str.auxil.τ_dd_layer,SPECTRA.IΛ_SIF,i  );    # transmittance of the layer without correction

        f_d_i = view(sun_geo.auxil.e_sifꜜ_layer    ,:,i  );            # downward emitted SIF from layer i
        f_u_i = view(sun_geo.auxil.e_sifꜛ_layer    ,:,i  );            # upward emitted SIF from layer i
        s_a_i = view(sun_geo.auxil.e_sifꜛ_layer_sum,:,i  );            # final upward SIF in the layer (sum of upward and transmitted downward SIF)
        s_d_i = view(sun_geo.auxil.e_sifꜜ_emit     ,:,i  );            # downward SIF from the layer
        s_u_i = view(sun_geo.auxil.e_sifꜛ_emit     ,:,i  );            # upward SIF from the layer
        s_u_j = view(sun_geo.auxil.e_sifꜛ_emit     ,:,i+1);            # upward SIF from the lower layer

        s_d_i .= (f_d_i .+ s_u_j .* r__) ./ (1 .- r__ .* r_j);
        s_u_i .= f_u_i .+ s_u_j .* t__ .+ s_d_i .* r_j .* t__;
        s_a_i .= f_u_i .+ s_d_i .* r_j .* t__;
    end;

    # 3. account for the SIF emission from up to bottom
    sun_geo.auxil.e_sifꜜ[:,1] .= 0;
    for i in 1:n_layer
        r_i = view(can_str.auxil.ρ_dd,SPECTRA.IΛ_SIF,i);    # reflectance of the layer (i) for SIF
        t_i = view(can_str.auxil.τ_dd,SPECTRA.IΛ_SIF,i);    # transmittance of the layer (i) for SIF

        s_d_i = view(sun_geo.auxil.e_sifꜜ_emit,:,i  );      # downward SIF from the layer
        s_u_i = view(sun_geo.auxil.e_sifꜛ_emit,:,i  );      # upward SIF from the layer
        a_d_i = view(sun_geo.auxil.e_sifꜜ     ,:,i  );
        a_d_j = view(sun_geo.auxil.e_sifꜜ     ,:,i+1);
        a_u_i = view(sun_geo.auxil.e_sifꜛ     ,:,i  );

        a_d_j .= a_d_i .* t_i .+ s_d_i;
        a_u_i .= a_d_i .* r_i .+ s_u_i;
    end;
    sun_geo.auxil.e_sifꜛ[:,end] .= view(sun_geo.auxil.e_sifꜜ,:,n_layer+1) .* view(can_str.auxil.ρ_dd,SPECTRA.IΛ_SIF,n_layer+1);
    sen_geo.auxil.sif_scattered .= view(sen_geo.auxil.dob_leaf,SPECTRA.IΛ_SIF,:) .* view(sun_geo.auxil.e_sifꜜ,:,1:n_layer) .+
                                   view(sen_geo.auxil.dof_leaf,SPECTRA.IΛ_SIF,:) .* view(sun_geo.auxil.e_sifꜛ,:,1:n_layer);

    # 4. compute SIF from the observer direction (CI is accounted for in the p_sensor and p_sun_sensor already, so do NOT use CI here)
    #    TODO: may have numerical issues because of due to the same issue with SIF conservation (might not, I do not know yet)
    vec_layer = spac.cache.cache_layer_1;
    vec_layer .= sen_geo.auxil.p_sun_sensor .* can_str.trait.δlai ./ FT(π);
    mul!(sen_geo.auxil.sif_obs_sunlit, sen_geo.auxil.sif_sunlit, vec_layer);

    vec_layer .= (sen_geo.auxil.p_sensor .- sen_geo.auxil.p_sun_sensor) .* can_str.trait.δlai ./ FT(π);
    mul!(sen_geo.auxil.sif_obs_shaded, sen_geo.auxil.sif_shaded, vec_layer);

    vec_layer .= sen_geo.auxil.p_sensor .* can_str.trait.δlai ./ FT(π);
    mul!(sen_geo.auxil.sif_obs_scattered, sen_geo.auxil.sif_scattered, vec_layer);

    sen_geo.auxil.sif_obs_soil .= view(sun_geo.auxil.e_sifꜛ,:,n_layer+1) .* sen_geo.auxil.p_sensor_soil ./ FT(π);

    sen_geo.auxil.sif_obs .= sen_geo.auxil.sif_obs_sunlit .+ sen_geo.auxil.sif_obs_shaded .+ sen_geo.auxil.sif_obs_scattered .+ sen_geo.auxil.sif_obs_soil;

    return nothing
end;


function lidf_weight end;

lidf_weight(mat_0::Matrix{FT}, p_incl::Vector{FT}, vec_azi::Vector{FT}) where {FT} = (
    mul!(vec_azi, mat_0', p_incl);

    # Note that because azimuth angle is evenly distributed, so we return the mean value here; otherwise, we will need to return p_azi' * vec_azi
    return mean(vec_azi)
);

lidf_weight(mat_prod::Matrix{FT}, mat_0::Matrix{FT}, mat_1::Matrix{FT}, p_incl::Vector{FT}, vec_azi::Vector{FT}) where {FT} = (
    mat_prod .= mat_0 .* mat_1;
    mul!(vec_azi, mat_prod', p_incl);

    # Note that because azimuth angle is evenly distributed, so we return the mean value here; otherwise, we will need to return p_azi' * vec_azi
    return mean(vec_azi)
);

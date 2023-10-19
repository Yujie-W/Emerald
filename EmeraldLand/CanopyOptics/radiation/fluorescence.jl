# This file contains functions to compute the SIF emission of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-14: add function fluorescence_spectrum! (run per sensor geometry)
#     2023-Oct-14: if LAI < = 0 or SZA > 89, set all fluxes to 0
#     2023-Oct-18: SIF excitation is rescaled to leaf partitioning (accounting stem)
#
#######################################################################################################################################################################################################
"""

    fluorescence_spectrum!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Compute the fluorescence spectrum of the canopy at the sensor direction, given
- `config` SPAC configuration
- `spac` SPAC

"""
function fluorescence_spectrum!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    if !config.ENABLE_SIF
        return nothing
    end;

    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    sun_geo = spac.canopy.sun_geometry;

    if sun_geo.state.sza > 89 || can_str.state.lai <= 0
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
    (; DIM_AZI, DIM_LAYER, SPECTRA, Φ_PHOTON) = config;

    # 0. compute chloroplast SIF emissions for different layers
    for i in 1:DIM_LAYER
        leaf = leaves[DIM_LAYER + 1 - i];
        a_leaf = view(leaf.bio.auxil.α_leaf,SPECTRA.IΛ_SIFE).* can_str.state.δlai[i];
        a_stem = (1 .- view(SPECTRA.ρ_STEM,SPECTRA.IΛ_SIFE)) .* can_str.state.δsai[i];
        f_leaf = a_leaf ./ (a_leaf .+ a_stem);

        # integrate the energy absorbed by chl (and car) in each wave length bins
        f_sife = view(leaf.bio.auxil.f_sife, SPECTRA.IΛ_SIFE);
        sun_geo.auxil._e_dif_sife .= view(sun_geo.auxil.e_net_dif,SPECTRA.IΛ_SIFE,i) .* f_leaf .* SPECTRA.ΔΛ_SIFE .* f_sife;
        sun_geo.auxil._e_dir_sife .= view(sun_geo.auxil.e_net_dir,SPECTRA.IΛ_SIFE,i) .* f_leaf .* SPECTRA.ΔΛ_SIFE .* f_sife;

        # convert the excitation radiation to fluorescence components
        if Φ_PHOTON
            photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_dif_sife);
            photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_dir_sife);
        end;

        # convert the excitation radiation to fluorescence components
        sun_geo.auxil._e_dif_sif .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(sun_geo.auxil._e_dif_sife);
        sun_geo.auxil._e_dir_sif .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(sun_geo.auxil._e_dir_sife);

        # add up the excitation radiation from direct and diffuse radiation for sunlit and shaded leaves
        sun_geo.auxil._e_dif_shaded .= sun_geo.auxil._e_dif_sif .* (1 - sun_geo.auxil.p_sunlit[i]);
        sun_geo.auxil._e_dif_sunlit .= sun_geo.auxil._e_dif_sif .* sun_geo.auxil.p_sunlit[i];
        sun_geo.auxil._e_dir_sunlit .= sun_geo.auxil._e_dir_sif;

        # convert the SIF back to energy unit if ϕ_photon is true
        if Φ_PHOTON
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dif_shaded);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dif_sunlit);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dir_sunlit);
        end;

        # add up the SIF from sunlit and shaded leaves for each layer through accounting for the SIF quantum yield
        ϕ_sunlit_dif = lidf_weight(leaf.flux.auxil.ϕ_f_sunlit, can_str.auxil.p_incl, sun_geo.auxil._vec_azi);
        ϕ_sunlit_dir = lidf_weight(sun_geo.auxil._mat_incl_azi, leaf.flux.auxil.ϕ_f_sunlit, sun_geo.auxil.fs_abs, sun_geo.auxil._vec_azi, can_str.auxil.p_incl);
        sun_geo.auxil.e_sif_chl[:,i] .= sun_geo.auxil._e_dif_shaded .* leaf.flux.auxil.ϕ_f_shaded .+ sun_geo.auxil._e_dif_sunlit .* ϕ_sunlit_dif .+ sun_geo.auxil._e_dir_sunlit .* ϕ_sunlit_dir;
    end;

    # 1. compute SIF emissions for different layers

    # function to weight matrices by inclination angles
    @inline local_lidf_weight(mat_0, mat_1) = (
        sun_geo.auxil._mat_incl_azi .= mat_0 .* mat_1;
        mul!(sun_geo.auxil._vec_azi, sun_geo.auxil._mat_incl_azi', can_str.auxil.p_incl);

        return mean(sun_geo.auxil._vec_azi)
    );
    _COS²_Θ_INCL_AZI = (cosd.(config.Θ_INCL) .^ 2) * ones(FT, 1, DIM_AZI);

    for i in 1:DIM_LAYER
        leaf = leaves[DIM_LAYER + 1 - i];
        a_leaf = view(leaf.bio.auxil.α_leaf,SPECTRA.IΛ_SIFE).* can_str.state.δlai[i];
        a_stem = (1 .- view(SPECTRA.ρ_STEM,SPECTRA.IΛ_SIFE)) .* can_str.state.δsai[i];
        f_leaf = a_leaf ./ (a_leaf .+ a_stem);

        # compute the energy used for SIF excitation
        sun_geo.auxil._e_dirꜜ_sife .= view(sun_geo.auxil.e_dirꜜ,SPECTRA.IΛ_SIFE,i) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜜ_sife .= view(sun_geo.auxil.e_difꜜ,SPECTRA.IΛ_SIFE,i) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜛ_sife .= view(sun_geo.auxil.e_difꜛ,SPECTRA.IΛ_SIFE,i) .* f_leaf .* SPECTRA.ΔΛ_SIFE;

        # convert the excitation radiation to photons if ϕ_photon is true
        if Φ_PHOTON
            photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_dirꜜ_sife);
            photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜜ_sife);
            photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜛ_sife);
        end;

        # convert the excitation radiation to fluorescence components
        mul!(sun_geo.auxil._e_dirꜜ_sif_mean, leaf.bio.auxil.mat_mean, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_dirꜜ_sif_diff, leaf.bio.auxil.mat_diff, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sif_mean, leaf.bio.auxil.mat_mean, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sif_diff, leaf.bio.auxil.mat_diff, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sif_mean, leaf.bio.auxil.mat_mean, sun_geo.auxil._e_difꜛ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sif_diff, leaf.bio.auxil.mat_diff, sun_geo.auxil._e_difꜛ_sife);

        # convert the SIF back to energy unit if ϕ_photon is true
        if Φ_PHOTON
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sif_mean);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sif_diff);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sif_mean);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sif_diff);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sif_mean);
            energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sif_diff);
        end;

        #
        #
        # TODO: refactor this part when fully understand what is happening here
        #
        #
        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        ϕ_sunlit = leaf.flux.auxil.ϕ_f_sunlit;
        ϕ_shaded = leaf.flux.auxil.ϕ_f_shaded;

        # compute the weights
        sl_1_ = local_lidf_weight(ϕ_sunlit, 1);
        sh_1_ = local_lidf_weight(ϕ_shaded, 1);
        sh_O_ = local_lidf_weight(ϕ_shaded, sen_geo.auxil.fo_abs);
        sl_O_ = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_abs);
        sl_S_ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_abs);
        sh_oθ = local_lidf_weight(ϕ_shaded, sen_geo.auxil.fo_cos²_incl);
        sl_oθ = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_cos²_incl);
        sl_sθ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_cos²_incl);
        sl_SO = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_fs_abs);
        sl_so = local_lidf_weight(ϕ_sunlit, sen_geo.auxil.fo_fs);
        sh_θ² = local_lidf_weight(ϕ_shaded, _COS²_Θ_INCL_AZI);
        sl_θ² = local_lidf_weight(ϕ_sunlit, _COS²_Θ_INCL_AZI);

        # upward and downward SIF from direct and diffuse radiation per leaf area
        sun_geo.auxil._sif_shadedꜛ .= sun_geo.auxil._e_difꜜ_sif_mean .* sh_1_ .+ sun_geo.auxil._e_difꜜ_sif_diff .* sh_θ² .+
                                      sun_geo.auxil._e_difꜛ_sif_mean .* sh_1_ .- sun_geo.auxil._e_difꜛ_sif_diff .* sh_θ²;
        sun_geo.auxil._sif_shadedꜜ .= sun_geo.auxil._e_difꜜ_sif_mean .* sh_1_ .- sun_geo.auxil._e_difꜜ_sif_diff .* sh_θ² .+
                                      sun_geo.auxil._e_difꜛ_sif_mean .* sh_1_ .+ sun_geo.auxil._e_difꜛ_sif_diff .* sh_θ²;
        sun_geo.auxil._sif_sunlitꜛ .= sun_geo.auxil._e_dirꜜ_sif_mean .* sl_S_ .+ sun_geo.auxil._e_dirꜜ_sif_diff .* sl_sθ .+
                                      sun_geo.auxil._e_difꜜ_sif_mean .* sl_1_ .+ sun_geo.auxil._e_difꜜ_sif_diff .* sl_θ² .+
                                      sun_geo.auxil._e_difꜛ_sif_mean .* sl_1_ .- sun_geo.auxil._e_difꜛ_sif_diff .* sl_θ²;
        sun_geo.auxil._sif_sunlitꜜ .= sun_geo.auxil._e_dirꜜ_sif_mean .* sl_S_ .- sun_geo.auxil._e_dirꜜ_sif_diff .* sl_sθ .+
                                      sun_geo.auxil._e_difꜜ_sif_mean .* sl_1_ .- sun_geo.auxil._e_difꜜ_sif_diff .* sl_θ² .+
                                      sun_geo.auxil._e_difꜛ_sif_mean .* sl_1_ .+ sun_geo.auxil._e_difꜛ_sif_diff .* sl_θ²;

        # update the SIF cache for the observer direction (compute it here to save time)
        sen_geo.auxil.sif_sunlit[:,i] .= sun_geo.auxil._e_dirꜜ_sif_mean .* sl_SO .+ sun_geo.auxil._e_dirꜜ_sif_diff .* sl_so .+
                                         sun_geo.auxil._e_difꜜ_sif_mean .* sl_O_ .+ sun_geo.auxil._e_dirꜜ_sif_diff .* sl_oθ .+
                                         sun_geo.auxil._e_difꜛ_sif_mean .* sl_O_ .- sun_geo.auxil._e_difꜛ_sif_diff .* sl_oθ;
        sen_geo.auxil.sif_shaded[:,i] .= sun_geo.auxil._e_difꜜ_sif_mean .* sh_O_ .+ sun_geo.auxil._e_dirꜜ_sif_diff .* sh_oθ .+
                                         sun_geo.auxil._e_difꜛ_sif_mean .* sh_O_ .- sun_geo.auxil._e_difꜛ_sif_diff .* sh_oθ;

        # total emitted SIF for upward and downward direction
        ilai = can_str.state.δlai[i] * can_str.auxil.ci;
        sun_geo.auxil.e_sifꜜ_layer[:,i] .= sun_geo.auxil._sif_sunlitꜜ .* ilai .* sun_geo.auxil.p_sunlit[i] .+ sun_geo.auxil._sif_shadedꜜ .* ilai .* (1 - sun_geo.auxil.p_sunlit[i]);
        sun_geo.auxil.e_sifꜛ_layer[:,i] .= sun_geo.auxil._sif_sunlitꜛ .* ilai .* sun_geo.auxil.p_sunlit[i] .+ sun_geo.auxil._sif_shadedꜛ .* ilai .* (1 - sun_geo.auxil.p_sunlit[i]);
    end;

    # 2. account for the SIF emission from bottom to up
    sun_geo.auxil.e_sifꜛ_emit[:,end] .= 0;
    for i in DIM_LAYER:-1:1
        r__ = view(can_str.auxil.ρ_dd_layer,SPECTRA.IΛ_SIF,i  );    # reflectance without correction
        r_j = view(can_str.auxil.ρ_dd      ,SPECTRA.IΛ_SIF,i+1);    # reflectance of the upper boundary (i) for SIF
        t_i = view(can_str.auxil.τ_dd      ,SPECTRA.IΛ_SIF,i  );    # transmittance of the layer (i) for SIF

        f_d_i = view(sun_geo.auxil.e_sifꜜ_layer,:,i  );            # downward emitted SIF from layer i
        f_u_i = view(sun_geo.auxil.e_sifꜛ_layer,:,i  );            # downward emitted SIF from layer i
        s_u_j = view(sun_geo.auxil.e_sifꜛ_emit ,:,i+1);            # upward SIF from the lower layer
        s_d_i = view(sun_geo.auxil.e_sifꜜ_emit ,:,i  );            # downward SIF from the layer
        s_u_i = view(sun_geo.auxil.e_sifꜛ_emit ,:,i  );            # upward SIF from the layer

        s_d_i .= (f_d_i .+ s_u_j .* r__) ./ (1 .- r__ .* r_j);
        s_u_i .= f_u_i .+ s_u_j .* t_i .+ s_d_i .* r_j .* t_i;
    end;

    # 3. account for the SIF emission from up to bottom
    sun_geo.auxil.e_sifꜜ[:,1] .= 0;
    for i in 1:DIM_LAYER
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
    sun_geo.auxil.e_sifꜛ[:,end] .= view(sun_geo.auxil.e_sifꜜ,:,DIM_LAYER+1) .* view(can_str.auxil.ρ_dd,SPECTRA.IΛ_SIF,DIM_LAYER+1);
    sen_geo.auxil.sif_scattered .= view(sen_geo.auxil.dob_leaf,SPECTRA.IΛ_SIF,:) .* view(sun_geo.auxil.e_sifꜜ,:,1:DIM_LAYER) .+
                                   view(sen_geo.auxil.dof_leaf,SPECTRA.IΛ_SIF,:) .* view(sun_geo.auxil.e_sifꜛ,:,1:DIM_LAYER);

    # 4. compute SIF from the observer direction
    vec_layer = ones(FT, DIM_LAYER);
    vec_layer .= sen_geo.auxil.p_sun_sensor .* can_str.state.δlai .* can_str.auxil.ci ./ FT(π);
    mul!(sen_geo.auxil.sif_obs_sunlit, sen_geo.auxil.sif_sunlit, vec_layer);

    vec_layer .= (sen_geo.auxil.p_sensor .- sen_geo.auxil.p_sun_sensor) .* can_str.state.δlai .* can_str.auxil.ci ./ FT(π);
    mul!(sen_geo.auxil.sif_obs_shaded, sen_geo.auxil.sif_shaded, vec_layer);

    vec_layer .= sen_geo.auxil.p_sensor .* can_str.state.δlai .* can_str.auxil.ci ./ FT(π);
    mul!(sen_geo.auxil.sif_obs_scattered, sen_geo.auxil.sif_scattered, vec_layer);

    sen_geo.auxil.sif_obs_soil .= view(sun_geo.auxil.e_sifꜛ,:,DIM_LAYER+1) .* sen_geo.auxil.p_sensor_soil ./ FT(π);

    sen_geo.auxil.sif_obs .= sen_geo.auxil.sif_obs_sunlit .+ sen_geo.auxil.sif_obs_shaded .+ sen_geo.auxil.sif_obs_scattered .+ sen_geo.auxil.sif_obs_soil;

    return nothing
end;


function lidf_weight end;

lidf_weight(mat_0::Matrix{FT}, p_incl::Vector{FT}, vec_azi::Vector{FT}) where {FT} = (
    mul!(vec_azi, mat_0', p_incl);

    # Note that because azimuth angle is evenly distributed, so we return the mean value here; otherwise, we will need to return p_azi' * vec_azi
    return mean(vec_azi)
);

lidf_weight(mat_prod::Matrix{FT}, mat_0::Matrix{FT}, mat_1::Matrix{FT}, vec_azi::Vector{FT}, p_incl::Vector{FT}) where {FT} = (
    mat_prod .= mat_0 .* mat_1;
    mul!(vec_azi, mat_prod', p_incl);

    # Note that because azimuth angle is evenly distributed, so we return the mean value here; otherwise, we will need to return p_azi' * vec_azi
    return mean(vec_azi)
);

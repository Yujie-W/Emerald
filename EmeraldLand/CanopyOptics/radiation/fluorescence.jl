# This file contains functions to compute the SIF emission of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-14: add function fluorescence_spectrum! (run per sensor geometry)
#     2023-Oct-14: if LAI < = 0 or SZA > 89, set all fluxes to 0
#
#######################################################################################################################################################################################################
"""

    fluorescence_spectrum!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the fluorescence spectrum of the canopy at the sensor direction, given
- `config` SPAC configuration
- `spac` SPAC

"""
function fluorescence_spectrum!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    if !config.ENABLE_SIF
        return nothing
    end;

    (; CANOPY, LEAVES) = spac;
    if spac.CANOPY.sun_geometry.state.sza > 89 || spac.CANOPY.structure.state.lai <= 0
        CANOPY.sun_geometry.auxil.e_sif_chl .= 0;
        CANOPY.sun_geometry.auxil.e_sifꜜ_layer .= 0;
        CANOPY.sun_geometry.auxil.e_sifꜛ_layer .= 0;
        CANOPY.sun_geometry.auxil.e_sifꜜ_emit .= 0;
        CANOPY.sun_geometry.auxil.e_sifꜛ_emit .= 0;
        CANOPY.sun_geometry.auxil.e_sifꜜ .= 0;
        CANOPY.sun_geometry.auxil.e_sifꜛ .= 0;
        CANOPY.sensor_geometry.auxil.sif_sunlit .= 0;
        CANOPY.sensor_geometry.auxil.sif_shaded .= 0;
        CANOPY.sensor_geometry.auxil.sif_scattered .= 0;
        CANOPY.sensor_geometry.auxil.sif_obs_sunlit .= 0;
        CANOPY.sensor_geometry.auxil.sif_obs_shaded .= 0;
        CANOPY.sensor_geometry.auxil.sif_obs_scattered .= 0;
        CANOPY.sensor_geometry.auxil.sif_obs_ssoil .= 0;
        CANOPY.sensor_geometry.auxil.sif_obs .= 0;

        return nothing
    end;

    # run the fluorescence simulations only if fluorescence feature is enabled
    (; DIM_AZI, DIM_LAYER, SPECTRA, Φ_PHOTON) = config;

    # 0. compute chloroplast SIF emissions for different layers
    for i in 1:DIM_LAYER
        j = DIM_LAYER + 1 - i;

        # integrate the energy absorbed by chl (and car) in each wave length bins
        f_sife = view(LEAVES[j].bio.auxil.f_sife, SPECTRA.IΛ_SIFE);
        CANOPY.sun_geometry.auxil._e_dif_sife .= view(CANOPY.sun_geometry.auxil.e_net_dif,SPECTRA.IΛ_SIFE,i) .* SPECTRA.ΔΛ_SIFE .* f_sife;
        CANOPY.sun_geometry.auxil._e_dir_sife .= view(CANOPY.sun_geometry.auxil.e_net_dir,SPECTRA.IΛ_SIFE,i) .* SPECTRA.ΔΛ_SIFE .* f_sife;

        # convert the excitation radiation to fluorescence components
        if Φ_PHOTON
            photon!(SPECTRA.Λ_SIFE, CANOPY.sun_geometry.auxil._e_dif_sife);
            photon!(SPECTRA.Λ_SIFE, CANOPY.sun_geometry.auxil._e_dir_sife);
        end;

        # convert the excitation radiation to fluorescence components
        CANOPY.sun_geometry.auxil._e_dif_sif .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(CANOPY.sun_geometry.auxil._e_dif_sife);
        CANOPY.sun_geometry.auxil._e_dir_sif .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(CANOPY.sun_geometry.auxil._e_dir_sife);

        # add up the excitation radiation from direct and diffuse radiation for sunlit and shaded leaves
        CANOPY.sun_geometry.auxil._e_dif_shaded .= CANOPY.sun_geometry.auxil._e_dif_sif .* (1 - CANOPY.sun_geometry.auxil.p_sunlit[i]);
        CANOPY.sun_geometry.auxil._e_dif_sunlit .= CANOPY.sun_geometry.auxil._e_dif_sif .* CANOPY.sun_geometry.auxil.p_sunlit[i];
        CANOPY.sun_geometry.auxil._e_dir_sunlit .= CANOPY.sun_geometry.auxil._e_dir_sif;

        # convert the SIF back to energy unit if ϕ_photon is true
        if Φ_PHOTON
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_dif_shaded);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_dif_sunlit);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_dir_sunlit);
        end;

        # add up the SIF from sunlit and shaded leaves for each layer through accounting for the SIF quantum yield
        ϕ_sunlit_dif = lidf_weight(LEAVES[j].flux.auxil.ϕ_f_sunlit,
                                   CANOPY.structure.state.p_incl,
                                   CANOPY.sun_geometry.auxil._vec_azi);
        ϕ_sunlit_dir = lidf_weight(CANOPY.sun_geometry.auxil._mat_incl_azi,
                                   LEAVES[j].flux.auxil.ϕ_f_sunlit,
                                   CANOPY.sun_geometry.auxil.fs_abs,
                                   CANOPY.sun_geometry.auxil._vec_azi,
                                   CANOPY.structure.state.p_incl);
        CANOPY.sun_geometry.auxil.e_sif_chl[:,i] .= CANOPY.sun_geometry.auxil._e_dif_shaded .* LEAVES[j].flux.auxil.ϕ_f_shaded .+
                                                    CANOPY.sun_geometry.auxil._e_dif_sunlit .* ϕ_sunlit_dif .+
                                                    CANOPY.sun_geometry.auxil._e_dir_sunlit .* ϕ_sunlit_dir;
    end;

    # 1. compute SIF emissions for different layers

    # function to weight matrices by inclination angles
    @inline local_lidf_weight(mat_0, mat_1) = (
        CANOPY.sun_geometry.auxil._mat_incl_azi .= mat_0 .* mat_1;
        mul!(CANOPY.sun_geometry.auxil._vec_azi, CANOPY.sun_geometry.auxil._mat_incl_azi', CANOPY.structure.state.p_incl);

        return mean(CANOPY.sun_geometry.auxil._vec_azi)
    );
    _COS²_Θ_INCL_AZI = (cosd.(config.Θ_INCL) .^ 2) * ones(FT, 1, DIM_AZI);

    for i in 1:DIM_LAYER
        j = DIM_LAYER + 1 - i;

        # compute the energy used for SIF excitation
        CANOPY.sun_geometry.auxil._e_dirꜜ_sife .= view(CANOPY.sun_geometry.auxil.e_dirꜜ,SPECTRA.IΛ_SIFE,i) .* SPECTRA.ΔΛ_SIFE;
        CANOPY.sun_geometry.auxil._e_difꜜ_sife .= view(CANOPY.sun_geometry.auxil.e_difꜜ,SPECTRA.IΛ_SIFE,i) .* SPECTRA.ΔΛ_SIFE;
        CANOPY.sun_geometry.auxil._e_difꜛ_sife .= view(CANOPY.sun_geometry.auxil.e_difꜛ,SPECTRA.IΛ_SIFE,i) .* SPECTRA.ΔΛ_SIFE;

        # convert the excitation radiation to photons if ϕ_photon is true
        if Φ_PHOTON
            photon!(SPECTRA.Λ_SIFE, CANOPY.sun_geometry.auxil._e_dirꜜ_sife);
            photon!(SPECTRA.Λ_SIFE, CANOPY.sun_geometry.auxil._e_difꜜ_sife);
            photon!(SPECTRA.Λ_SIFE, CANOPY.sun_geometry.auxil._e_difꜛ_sife);
        end;

        # convert the excitation radiation to fluorescence components
        mul!(CANOPY.sun_geometry.auxil._e_dirꜜ_sif_mean, LEAVES[j].bio.auxil.mat_mean, CANOPY.sun_geometry.auxil._e_dirꜜ_sife);
        mul!(CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff, LEAVES[j].bio.auxil.mat_diff, CANOPY.sun_geometry.auxil._e_dirꜜ_sife);
        mul!(CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean, LEAVES[j].bio.auxil.mat_mean, CANOPY.sun_geometry.auxil._e_difꜜ_sife);
        mul!(CANOPY.sun_geometry.auxil._e_difꜜ_sif_diff, LEAVES[j].bio.auxil.mat_diff, CANOPY.sun_geometry.auxil._e_difꜜ_sife);
        mul!(CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean, LEAVES[j].bio.auxil.mat_mean, CANOPY.sun_geometry.auxil._e_difꜛ_sife);
        mul!(CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff, LEAVES[j].bio.auxil.mat_diff, CANOPY.sun_geometry.auxil._e_difꜛ_sife);

        # convert the SIF back to energy unit if ϕ_photon is true
        if Φ_PHOTON
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_dirꜜ_sif_mean);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_difꜜ_sif_diff);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean);
            energy!(SPECTRA.Λ_SIF, CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff);
        end;

        #
        #
        # TODO: refactor this part when fully understand what is happening here
        #
        #
        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        ϕ_sunlit = LEAVES[j].flux.auxil.ϕ_f_sunlit;
        ϕ_shaded = LEAVES[j].flux.auxil.ϕ_f_shaded;

        # compute the weights
        sl_1_ = local_lidf_weight(ϕ_sunlit, 1);
        sh_1_ = local_lidf_weight(ϕ_shaded, 1);
        sh_O_ = local_lidf_weight(ϕ_shaded, CANOPY.sensor_geometry.auxil.fo_abs);
        sl_O_ = local_lidf_weight(ϕ_sunlit, CANOPY.sensor_geometry.auxil.fo_abs);
        sl_S_ = local_lidf_weight(ϕ_sunlit, CANOPY.sun_geometry.auxil.fs_abs);
        sh_oθ = local_lidf_weight(ϕ_shaded, CANOPY.sensor_geometry.auxil.fo_cos²_incl);
        sl_oθ = local_lidf_weight(ϕ_sunlit, CANOPY.sensor_geometry.auxil.fo_cos²_incl);
        sl_sθ = local_lidf_weight(ϕ_sunlit, CANOPY.sun_geometry.auxil.fs_cos²_incl);
        sl_SO = local_lidf_weight(ϕ_sunlit, CANOPY.sensor_geometry.auxil.fo_fs_abs);
        sl_so = local_lidf_weight(ϕ_sunlit, CANOPY.sensor_geometry.auxil.fo_fs);
        sh_θ² = local_lidf_weight(ϕ_shaded, _COS²_Θ_INCL_AZI);
        sl_θ² = local_lidf_weight(ϕ_sunlit, _COS²_Θ_INCL_AZI);

        # upward and downward SIF from direct and diffuse radiation per leaf area
        CANOPY.sun_geometry.auxil._sif_shadedꜛ .= CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean .* sh_1_ .+ CANOPY.sun_geometry.auxil._e_difꜜ_sif_diff .* sh_θ² .+
                                                  CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean .* sh_1_ .- CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff .* sh_θ²;
        CANOPY.sun_geometry.auxil._sif_shadedꜜ .= CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean .* sh_1_ .- CANOPY.sun_geometry.auxil._e_difꜜ_sif_diff .* sh_θ² .+
                                                  CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean .* sh_1_ .+ CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff .* sh_θ²;
        CANOPY.sun_geometry.auxil._sif_sunlitꜛ .= CANOPY.sun_geometry.auxil._e_dirꜜ_sif_mean .* sl_S_ .+ CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff .* sl_sθ .+
                                                  CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean .* sl_1_ .+ CANOPY.sun_geometry.auxil._e_difꜜ_sif_diff .* sl_θ² .+
                                                  CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean .* sl_1_ .- CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff .* sl_θ²;
        CANOPY.sun_geometry.auxil._sif_sunlitꜜ .= CANOPY.sun_geometry.auxil._e_dirꜜ_sif_mean .* sl_S_ .- CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff .* sl_sθ .+
                                                  CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean .* sl_1_ .- CANOPY.sun_geometry.auxil._e_difꜜ_sif_diff .* sl_θ² .+
                                                  CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean .* sl_1_ .+ CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff .* sl_θ²;

        # update the SIF cache for the observer direction (compute it here to save time)
        CANOPY.sensor_geometry.auxil.sif_sunlit[:,i] .= CANOPY.sun_geometry.auxil._e_dirꜜ_sif_mean .* sl_SO .+ CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff .* sl_so .+
                                                        CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean .* sl_O_ .+ CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff .* sl_oθ .+
                                                        CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean .* sl_O_ .- CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff .* sl_oθ;
        CANOPY.sensor_geometry.auxil.sif_shaded[:,i] .= CANOPY.sun_geometry.auxil._e_difꜜ_sif_mean .* sh_O_ .+ CANOPY.sun_geometry.auxil._e_dirꜜ_sif_diff .* sh_oθ .+
                                                        CANOPY.sun_geometry.auxil._e_difꜛ_sif_mean .* sh_O_ .- CANOPY.sun_geometry.auxil._e_difꜛ_sif_diff .* sh_oθ;

        # total emitted SIF for upward and downward direction
        ilai = CANOPY.structure.state.δlai[i] * CANOPY.structure.auxil.ci;
        CANOPY.sun_geometry.auxil.e_sifꜜ_layer[:,i] .= CANOPY.sun_geometry.auxil._sif_sunlitꜜ .* ilai .* CANOPY.sun_geometry.auxil.p_sunlit[i] .+
                                                       CANOPY.sun_geometry.auxil._sif_shadedꜜ .* ilai .* (1 - CANOPY.sun_geometry.auxil.p_sunlit[i]);
        CANOPY.sun_geometry.auxil.e_sifꜛ_layer[:,i] .= CANOPY.sun_geometry.auxil._sif_sunlitꜛ .* ilai .* CANOPY.sun_geometry.auxil.p_sunlit[i] .+
                                                       CANOPY.sun_geometry.auxil._sif_shadedꜛ .* ilai .* (1 - CANOPY.sun_geometry.auxil.p_sunlit[i]);
    end;

    # 2. account for the SIF emission from bottom to up
    CANOPY.sun_geometry.auxil.e_sifꜛ_emit[:,end] .= 0;
    for i in DIM_LAYER:-1:1
        r__ = view(CANOPY.structure.auxil.ρ_dd_layer,SPECTRA.IΛ_SIF,i  );   # reflectance without correction
        r_j = view(CANOPY.structure.auxil.ρ_dd      ,SPECTRA.IΛ_SIF,i+1);   # reflectance of the upper boundary (i) for SIF
        t_i = view(CANOPY.structure.auxil.τ_dd      ,SPECTRA.IΛ_SIF,i  );   # transmittance of the layer (i) for SIF

        f_d_i = view(CANOPY.sun_geometry.auxil.e_sifꜜ_layer,:,i  );         # downward emitted SIF from layer i
        f_u_i = view(CANOPY.sun_geometry.auxil.e_sifꜛ_layer,:,i  );         # downward emitted SIF from layer i
        s_u_j = view(CANOPY.sun_geometry.auxil.e_sifꜛ_emit ,:,i+1);         # upward SIF from the lower layer
        s_d_i = view(CANOPY.sun_geometry.auxil.e_sifꜜ_emit ,:,i  );         # downward SIF from the layer
        s_u_i = view(CANOPY.sun_geometry.auxil.e_sifꜛ_emit ,:,i  );         # upward SIF from the layer

        s_d_i .= (f_d_i .+ s_u_j .* r__) ./ (1 .- r__ .* r_j);
        s_u_i .= f_u_i .+ s_u_j .* t_i .+ s_d_i .* r_j .* t_i;
    end;

    # 3. account for the SIF emission from up to bottom
    CANOPY.sun_geometry.auxil.e_sifꜜ[:,1] .= 0;
    for i in 1:DIM_LAYER
        r_i = view(CANOPY.structure.auxil.ρ_dd,SPECTRA.IΛ_SIF,i);   # reflectance of the layer (i) for SIF
        t_i = view(CANOPY.structure.auxil.τ_dd,SPECTRA.IΛ_SIF,i);   # transmittance of the layer (i) for SIF

        s_d_i = view(CANOPY.sun_geometry.auxil.e_sifꜜ_emit,:,i  );  # downward SIF from the layer
        s_u_i = view(CANOPY.sun_geometry.auxil.e_sifꜛ_emit,:,i  );  # upward SIF from the layer
        a_d_i = view(CANOPY.sun_geometry.auxil.e_sifꜜ     ,:,i  );
        a_d_j = view(CANOPY.sun_geometry.auxil.e_sifꜜ     ,:,i+1);
        a_u_i = view(CANOPY.sun_geometry.auxil.e_sifꜛ     ,:,i  );

        a_d_j .= a_d_i .* t_i .+ s_d_i;
        a_u_i .= a_d_i .* r_i .+ s_u_i;
    end;
    CANOPY.sun_geometry.auxil.e_sifꜛ[:,end] .= view(CANOPY.sun_geometry.auxil.e_sifꜜ,:,DIM_LAYER+1) .* view(CANOPY.structure.auxil.ρ_dd,SPECTRA.IΛ_SIF,DIM_LAYER+1);
    CANOPY.sensor_geometry.auxil.sif_scattered .= view(CANOPY.sensor_geometry.auxil.dob_leaf,SPECTRA.IΛ_SIF,:) .* view(CANOPY.sun_geometry.auxil.e_sifꜜ,:,1:DIM_LAYER) .+
                                                  view(CANOPY.sensor_geometry.auxil.dof_leaf,SPECTRA.IΛ_SIF,:) .* view(CANOPY.sun_geometry.auxil.e_sifꜛ,:,1:DIM_LAYER);

    # 4. compute SIF from the observer direction
    vec_layer = ones(FT, DIM_LAYER);
    vec_layer .= CANOPY.sensor_geometry.auxil.p_sun_sensor .* CANOPY.structure.state.δlai .* CANOPY.structure.auxil.ci ./ FT(π);
    mul!(CANOPY.sensor_geometry.auxil.sif_obs_sunlit, CANOPY.sensor_geometry.auxil.sif_sunlit, vec_layer);

    vec_layer .= (CANOPY.sensor_geometry.auxil.p_sensor .- CANOPY.sensor_geometry.auxil.p_sun_sensor) .* CANOPY.structure.state.δlai .* CANOPY.structure.auxil.ci ./ FT(π);
    mul!(CANOPY.sensor_geometry.auxil.sif_obs_shaded, CANOPY.sensor_geometry.auxil.sif_shaded, vec_layer);

    vec_layer .= CANOPY.sensor_geometry.auxil.p_sensor .* CANOPY.structure.state.δlai .* CANOPY.structure.auxil.ci ./ FT(π);
    mul!(CANOPY.sensor_geometry.auxil.sif_obs_scattered, CANOPY.sensor_geometry.auxil.sif_scattered, vec_layer);

    CANOPY.sensor_geometry.auxil.sif_obs_soil .= view(CANOPY.sun_geometry.auxil.e_sifꜛ,:,DIM_LAYER+1) .* CANOPY.sensor_geometry.auxil.p_sensor_soil ./ FT(π);

    CANOPY.sensor_geometry.auxil.sif_obs .= CANOPY.sensor_geometry.auxil.sif_obs_sunlit .+
                                            CANOPY.sensor_geometry.auxil.sif_obs_shaded .+
                                            CANOPY.sensor_geometry.auxil.sif_obs_scattered .+
                                            CANOPY.sensor_geometry.auxil.sif_obs_soil;

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

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-10: migrate the function from CanopyLayers
#     2022-Jun-14: convert energy and photon back and forth if using photon mode
#     2022-Jun-29: use ϕ_f in Leaves2D
#     2022-Jun-29: add method for SPAC
#     2023-Mar-11: compute fluorescence only if solar zenith angle < 89
#     2023-Mar-11: add code to account for the case of LAI == 0
#     2023-May-19: use δlai per canopy layer
#     2023-Sep-11: compute the SIF at chloroplast level at the same time
#     2023-Sep-11: redo the calculation of SIF for different layers using the SIF excitation radiation directly
#     2023-Sep-19: recompute the SIF at chlorophyll level using the net absorbed radiation
# Bug fixes
#     2023-Mar-16: ddb ddf to dob and dof for observed SIF
#     2023-Sep-11: set ilai to lai rather than lai*ci
#
#######################################################################################################################################################################################################
"""

    canopy_fluorescence!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Updates canopy fluorescence, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` type SPAC

"""
function canopy_fluorescence! end

canopy_fluorescence!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; ANGLES, CANOPY, LEAVES) = spac;

    if (ANGLES.sza < 89)
        canopy_fluorescence!(config, CANOPY, LEAVES);
    else
        CANOPY.RADIATION.sif_obs .= 0;
    end;

    return nothing
);

canopy_fluorescence!(config::SPACConfiguration{FT}, can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}) where {FT} = (
    (; DIM_LAYER, SPECTRA, Φ_PHOTON, _COS²_Θ_INCL_AZI) = config;
    (; OPTICS, P_INCL, RADIATION) = can;

    if can.lai == 0
        RADIATION.sif_obs .= 0;
    end;

    # function to weight matrices by inclination angles
    @inline lidf_weight(mat_0, mat_1) = (
        OPTICS._tmp_mat_incl_azi_1 .= mat_0 .* mat_1;
        mul!(OPTICS._tmp_vec_azi, OPTICS._tmp_mat_incl_azi_1', P_INCL);

        return mean(OPTICS._tmp_vec_azi)
    );

    # 0. compute chloroplast SIF emissions for different layers
    for _i in 1:DIM_LAYER
        _k_chl = view(leaves[_i].NS.bio.auxil.f_sife, SPECTRA.IΛ_SIFE);

        # integrate the energy absorbed by chl (and car) in each wave length bins
        OPTICS._tmp_vec_sife_1 .= view(RADIATION.e_net_diffuse,SPECTRA.IΛ_SIFE,_i) .* SPECTRA.ΔΛ_SIFE .* _k_chl;
        OPTICS._tmp_vec_sife_2 .= view(RADIATION.e_net_direct ,SPECTRA.IΛ_SIFE,_i) .* SPECTRA.ΔΛ_SIFE .* _k_chl;

        # determine which ones to use depending on ϕ_photon
        if Φ_PHOTON
            photon!(SPECTRA.Λ_SIFE, OPTICS._tmp_vec_sife_1);
            photon!(SPECTRA.Λ_SIFE, OPTICS._tmp_vec_sife_2);
        end;

        _e_dir, _e_dif = OPTICS._tmp_vec_sife_1, OPTICS._tmp_vec_sife_2;

        # convert the excitation radiation to fluorescence components
        OPTICS._tmp_vec_sif_1 .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(_e_dir);
        OPTICS._tmp_vec_sif_2 .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(_e_dif) * OPTICS.p_sunlit[_i];
        OPTICS._tmp_vec_sif_3 .= view(SPECTRA.Φ_PS,SPECTRA.IΛ_SIF) .* sum(_e_dif) * (1 - OPTICS.p_sunlit[_i]);

        # convert the SIF back to energy unit if ϕ_photon is true
        if Φ_PHOTON
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_1);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_2);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_3);
        end;

        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        _ϕ_sunlit = leaves[_i].ϕ_f_sunlit;
        _ϕ_shaded = leaves[_i].ϕ_f_shaded;

        # compute the weights
        _sh_1_ = lidf_weight(_ϕ_shaded, 1);
        _sl_1_ = lidf_weight(_ϕ_sunlit, 1);
        _sl_S_ = lidf_weight(_ϕ_sunlit, OPTICS._abs_fs);

        # upward and downward SIF from direct and diffuse radiation per leaf area
        RADIATION._s_sunlit_up .= OPTICS._tmp_vec_sif_2 .* _sl_1_ .+ OPTICS._tmp_vec_sif_1 .* _sl_S_;
        RADIATION._s_shaded_up .= OPTICS._tmp_vec_sif_3 .* _sh_1_;

        # total emitted SIF for upward and downward direction
        # set _ilai to LAI but not LAI * CI as this is the total emitted SIF
        RADIATION.s_layer_down_chl[:,_i] .= 0;
        RADIATION.s_layer_up_chl[:,_i] .= RADIATION._s_sunlit_up .+ RADIATION._s_shaded_up;
    end;

    # 1. compute SIF emissions for different layers
    for _i in 1:DIM_LAYER
        OPTICS._mat⁺ .= (leaves[_i].NS.bio.auxil.mat_b .+ leaves[_i].NS.bio.auxil.mat_f) ./ 2;
        OPTICS._mat⁻ .= (leaves[_i].NS.bio.auxil.mat_b .- leaves[_i].NS.bio.auxil.mat_f) ./ 2;

        # integrate the energy in each wave length bins
        OPTICS._tmp_vec_sife_1 .= view(RADIATION.e_direct      ,SPECTRA.IΛ_SIFE,1 ) .* SPECTRA.ΔΛ_SIFE;
        OPTICS._tmp_vec_sife_2 .= view(RADIATION.e_diffuse_down,SPECTRA.IΛ_SIFE,_i) .* SPECTRA.ΔΛ_SIFE;
        OPTICS._tmp_vec_sife_3 .= view(RADIATION.e_diffuse_up  ,SPECTRA.IΛ_SIFE,_i) .* SPECTRA.ΔΛ_SIFE;

        # determine which ones to use depending on ϕ_photon
        if Φ_PHOTON
            photon!(SPECTRA.Λ_SIFE, OPTICS._tmp_vec_sife_1);
            photon!(SPECTRA.Λ_SIFE, OPTICS._tmp_vec_sife_2);
            photon!(SPECTRA.Λ_SIFE, OPTICS._tmp_vec_sife_3);
        end;

        _e_dir, _e_dif_down, _e_dif_up = OPTICS._tmp_vec_sife_1, OPTICS._tmp_vec_sife_2, OPTICS._tmp_vec_sife_3;

        # convert the excitation radiation to fluorescence components
        mul!(OPTICS._tmp_vec_sif_1, OPTICS._mat⁺, _e_dir);          # SIF component from direct light (before scaling)
        mul!(OPTICS._tmp_vec_sif_2, OPTICS._mat⁻, _e_dir);          # SIF component from direct light (before scaling)
        mul!(OPTICS._tmp_vec_sif_3, OPTICS._mat⁺, _e_dif_down);     # SIF component from downward diffuse light for backward (before scaling)
        mul!(OPTICS._tmp_vec_sif_4, OPTICS._mat⁻, _e_dif_down);     # SIF component from downward diffuse light for backward (before scaling)
        mul!(OPTICS._tmp_vec_sif_5, OPTICS._mat⁺, _e_dif_up);       # SIF component from upward diffuse light for forward (before scaling)
        mul!(OPTICS._tmp_vec_sif_6, OPTICS._mat⁻, _e_dif_up);       # SIF component from upward diffuse light for forward (before scaling)

        # convert the SIF back to energy unit if ϕ_photon is true
        if Φ_PHOTON
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_1);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_2);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_3);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_4);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_5);
            energy!(SPECTRA.Λ_SIF, OPTICS._tmp_vec_sif_6);
        end;

        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        _ϕ_sunlit = leaves[_i].ϕ_f_sunlit;
        _ϕ_shaded = leaves[_i].ϕ_f_shaded;

        # compute the weights
        _sh_1_ = lidf_weight(_ϕ_shaded, 1);
        _sl_1_ = lidf_weight(_ϕ_sunlit, 1);
        _sh_O_ = lidf_weight(_ϕ_shaded, OPTICS._abs_fo);
        _sl_O_ = lidf_weight(_ϕ_sunlit, OPTICS._abs_fo);
        _sl_S_ = lidf_weight(_ϕ_sunlit, OPTICS._abs_fs);
        _sh_oθ = lidf_weight(_ϕ_shaded, OPTICS._fo_cos_θ_incl);
        _sl_oθ = lidf_weight(_ϕ_sunlit, OPTICS._fo_cos_θ_incl);
        _sl_sθ = lidf_weight(_ϕ_sunlit, OPTICS._fs_cos_θ_incl);
        _sl_SO = lidf_weight(_ϕ_sunlit, OPTICS._abs_fs_fo);
        _sl_so = lidf_weight(_ϕ_sunlit, OPTICS._fs_fo);
        _sh_θ² = lidf_weight(_ϕ_shaded, _COS²_Θ_INCL_AZI);
        _sl_θ² = lidf_weight(_ϕ_sunlit, _COS²_Θ_INCL_AZI);

        # upward and downward SIF from direct and diffuse radiation per leaf area
        RADIATION._s_shaded_up   .= OPTICS._tmp_vec_sif_3 .* _sh_1_ .+ OPTICS._tmp_vec_sif_4 .* _sh_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sh_1_ .- OPTICS._tmp_vec_sif_6 .* _sh_θ²;
        RADIATION._s_shaded_down .= OPTICS._tmp_vec_sif_3 .* _sh_1_ .- OPTICS._tmp_vec_sif_4 .* _sh_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sh_1_ .+ OPTICS._tmp_vec_sif_6 .* _sh_θ²;
        RADIATION._s_sunlit_up   .= OPTICS._tmp_vec_sif_1 .* _sl_S_ .+ OPTICS._tmp_vec_sif_2 .* _sl_sθ .+
                                    OPTICS._tmp_vec_sif_3 .* _sl_1_ .+ OPTICS._tmp_vec_sif_4 .* _sl_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sl_1_ .- OPTICS._tmp_vec_sif_6 .* _sl_θ²;
        RADIATION._s_sunlit_down .= OPTICS._tmp_vec_sif_1 .* _sl_S_ .- OPTICS._tmp_vec_sif_2 .* _sl_sθ .+
                                    OPTICS._tmp_vec_sif_3 .* _sl_1_ .- OPTICS._tmp_vec_sif_4 .* _sl_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sl_1_ .+ OPTICS._tmp_vec_sif_6 .* _sl_θ²;

        # update the SIF cache for the observer direction (compute it here to save time)
        RADIATION._sif_obs_sunlit[:,_i] .= OPTICS._tmp_vec_sif_1 .* _sl_SO .+ OPTICS._tmp_vec_sif_2 .* _sl_so .+
                                           OPTICS._tmp_vec_sif_3 .* _sl_O_ .+ OPTICS._tmp_vec_sif_2 .* _sl_oθ .+
                                           OPTICS._tmp_vec_sif_5 .* _sl_O_ .- OPTICS._tmp_vec_sif_6 .* _sl_oθ;
        RADIATION._sif_obs_shaded[:,_i] .= OPTICS._tmp_vec_sif_3 .* _sh_O_ .+ OPTICS._tmp_vec_sif_2 .* _sh_oθ .+
                                           OPTICS._tmp_vec_sif_5 .* _sh_O_ .- OPTICS._tmp_vec_sif_6 .* _sh_oθ;

        # total emitted SIF for upward and downward direction
        # set _ilai to LAI but not LAI * CI as this is the total emitted SIF
        _ilai = can.δlai[_i];
        RADIATION.s_layer_down[:,_i] .= _ilai .* OPTICS.p_sunlit[_i] .* RADIATION._s_sunlit_down .+ _ilai .* (1 - OPTICS.p_sunlit[_i]) .* RADIATION._s_shaded_down;
        RADIATION.s_layer_up[:,_i]   .= _ilai .* OPTICS.p_sunlit[_i] .* RADIATION._s_sunlit_up   .+ _ilai .* (1 - OPTICS.p_sunlit[_i]) .* RADIATION._s_shaded_up;
    end;

    # 2. account for the SIF emission from bottom to up
    RADIATION._s_emit_up[:,end] .= 0;

    for _i in DIM_LAYER:-1:1
        _r__ = view(OPTICS._ρ_dd,SPECTRA.IΛ_SIF,_i  );  # reflectance without correction
        _r_j = view(OPTICS.ρ_dd ,SPECTRA.IΛ_SIF,_i+1);  # reflectance of the upper boundary (i) for SIF
        _t_i = view(OPTICS.τ_dd ,SPECTRA.IΛ_SIF,_i  );  # transmittance of the layer (i) for SIF

        _s_d_i = view(RADIATION._s_emit_down,:,_i  );   # downward SIF from the layer
        _s_u_i = view(RADIATION._s_emit_up  ,:,_i  );   # upward SIF from the layer
        _s_u_j = view(RADIATION._s_emit_up  ,:,_i+1);   # upward SIF from the lower layer
        _f_d_i = view(RADIATION.s_layer_down,:,_i  );   # downward emitted SIF from layer i
        _f_u_i = view(RADIATION.s_layer_up  ,:,_i  );   # downward emitted SIF from layer i

        OPTICS._tmp_vec_sif_1 .= 1 .- _r__ .* _r_j;

        _s_d_i .= (_f_d_i .+ _s_u_j .* _r__) ./ OPTICS._tmp_vec_sif_1;
        _s_u_i .= _f_u_i .+ _s_u_j .* _t_i .+ _s_d_i .* _r_j .* _t_i;
    end;

    # 3. account for the SIF emission from up to bottom
    RADIATION.sif_down[:,1] .= 0;

    for _i in 1:DIM_LAYER
        _r_i = view(OPTICS.ρ_dd,SPECTRA.IΛ_SIF,_i);   # reflectance of the layer (i) for SIF
        _t_i = view(OPTICS.τ_dd,SPECTRA.IΛ_SIF,_i);   # transmittance of the layer (i) for SIF

        _s_d_i = view(RADIATION._s_emit_down,:,_i);     # downward SIF from the layer
        _s_u_i = view(RADIATION._s_emit_up  ,:,_i);     # upward SIF from the layer

        _a_d_i = view(RADIATION.sif_down,:,_i  );
        _a_d_j = view(RADIATION.sif_down,:,_i+1);
        _a_u_i = view(RADIATION.sif_up  ,:,_i  );

        _a_d_j .= _a_d_i .* _t_i .+ _s_d_i;
        _a_u_i .= _a_d_i .* _r_i .+ _s_u_i;
    end;

    RADIATION.sif_up[:,end] .= view(RADIATION.sif_down,:,DIM_LAYER+1) .* view(OPTICS.ρ_dd,SPECTRA.IΛ_SIF,DIM_LAYER+1) .+ view(RADIATION._s_emit_up,:,DIM_LAYER+1);

    # 4. compute SIF from the observer direction
    OPTICS._tmp_vec_layer .= (view(OPTICS.pso,1:DIM_LAYER) .+ view(OPTICS.pso,2:DIM_LAYER+1)) ./ 2 .* can.δlai .* can.ci ./ FT(pi);
    mul!(RADIATION.sif_obs_sunlit, RADIATION._sif_obs_sunlit, OPTICS._tmp_vec_layer);

    OPTICS._tmp_vec_layer .= (view(OPTICS.po,1:DIM_LAYER) .+ view(OPTICS.po,2:DIM_LAYER+1) .- view(OPTICS.pso,1:DIM_LAYER) .- view(OPTICS.pso,2:DIM_LAYER+1)) ./ 2 .* can.δlai .* can.ci ./ FT(pi);
    mul!(RADIATION.sif_obs_shaded, RADIATION._sif_obs_shaded, OPTICS._tmp_vec_layer);

    RADIATION._sif_obs_scatter .= view(OPTICS.σ_dob,SPECTRA.IΛ_SIF,:) .* view(RADIATION.sif_down,:,1:DIM_LAYER) .+ view(OPTICS.σ_dof,SPECTRA.IΛ_SIF,:) .* view(RADIATION.sif_up,:,1:DIM_LAYER);
    OPTICS._tmp_vec_layer .= (view(OPTICS.po,1:DIM_LAYER) .+ view(OPTICS.po,2:DIM_LAYER+1)) ./ 2 .* can.δlai .* can.ci ./ FT(pi);
    mul!(RADIATION.sif_obs_scatter, RADIATION._sif_obs_scatter, OPTICS._tmp_vec_layer);

    RADIATION.sif_obs_soil .= view(RADIATION.sif_up,:,DIM_LAYER+1) .* OPTICS.po[end] ./ FT(pi);

    RADIATION.sif_obs .= RADIATION.sif_obs_sunlit .+ RADIATION.sif_obs_shaded .+ RADIATION.sif_obs_scatter .+ RADIATION.sif_obs_soil;

    return nothing
);
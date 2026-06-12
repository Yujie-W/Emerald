# This file contains functions to compute the sensor geometry of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-22: add function sensor_geometry_aux! to update the state-dependent auxiliary variables for sensor geometry (to call in step_remote_sensing!)
#     2024-Mar-01: compute the layer shortwave scattering coefficients based on the new theory
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: redesign the pso equation to account for the dngular dependence of clumping index
#     2024-Sep-09: account for diffuse CI impact on ko
# Bug fixes
#     2024-Mar-06: ci impact on fraction from viewer direction
#     2024-Sep-07: fix the calculation of sensa.ko_incl (was using [i] = 2 / FT(π) / cosd(FT(pi)) ...)
#
#######################################################################################################################################################################################################
"""

    sensor_geometry_aux!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sensor_geometry_aux! end;

sensor_geometry_aux!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} =
    sensor_geometry_aux!(config, spac.canopy, min(FT(0.5), spac.plant.leaves[1].bio.trait.width / (spac.plant.zs[2] - spac.plant.zs[1])));

sensor_geometry_aux!(config::SPACConfig{FT}, can::MultiLayerCanopy{FT}, lw2ch::FT) where {FT} =
    sensor_geometry_aux!(
            config,
            config.METHODS.CANOPY_RT_METHOD,
            can.structure.trait,
            can.structure.auxil,
            can.sun_geometry.state,
            can.sun_geometry.auxil,
            can.sensor_geometry.state,
            can.sensor_geometry.auxil,
            lw2ch);

sensor_geometry_aux!(
            config::SPACConfig{FT},
            ::CanopyRTEmerald,
            canst::CanopyStructureTrait{FT},
            cansa::CanopyStructureAuxil{FT},
            sunst::SunGeometryState{FT},
            sunsa::SunGeometryAuxil{FT},
            senst::SensorGeometryState{FT},
            sensa::SensorGeometryAuxil{FT},
            lw2ch::FT) where {FT} = (
    # if none of REF or SIF is enabled, or sza > 89, or LAI+SAI <= 0, do nothing
    if (!config.FEATURES.ENABLE_REF && !config.FEATURES.ENABLE_SIF) || sunst.sza > 89 || (canst.lai <= 0 && canst.sai <= 0)
        return nothing
    end;

    (; Θ_AZI, Θ_INCL) = config.DIMENSIONS;

    # compute clumping index from sensor zenith angle
    sensa.ci_sensor = canst.ci.ci_0 * (1 - canst.ci.ci_1 * cosd(senst.vza));

    # extinction coefficients for the solar radiation
    vza = senst.vza;
    sza = sunst.sza;
    raa = senst.vaa - sunst.saa;
    for i in eachindex(Θ_INCL)
        # observer angle
        Co = cosd(Θ_INCL[i]) * cosd(vza);
        So = sind(Θ_INCL[i]) * sind(vza);
        cosβo = abs(So) <= 1e-6 ? FT(1) : -Co/So;
        if abs(cosβo) < 1
            βo = acos(cosβo);
            Do = So;
        elseif vza < 90
            βo = FT(π);
            Do = Co;
        else
            βo = 0;
            Do = -Co;
        end;
        sensa.Co_incl[i] = Co;
        sensa.So_incl[i] = So;
        sensa.βo_incl[i] = βo;
        sensa.ko_incl[i] = 2 / FT(π) / cosd(vza) * (Co * (βo - FT(π)/2) + So * sin(βo));

        # compute the scattering coefficients
        Cs = sunsa.Cs_incl[i];
        Ss = sunsa.Ss_incl[i];
        βs = sunsa.βs_incl[i];
        Ds = (abs(cos(βs)) < 1 ? Ss : Cs);

        # 1 compute the Δ and β angles
        Δ₁ = abs(βs - βo);
        Δ₂ = FT(π) - abs(βs + βo - FT(π));

        ψ = deg2rad( abs(raa - 360*round(raa/360)) );
        if ψ <= Δ₁
            β₁,β₂,β₃ = ψ,Δ₁,Δ₂;
        elseif Δ₁ < ψ < Δ₂
            β₁,β₂,β₃ = Δ₁,ψ,Δ₂;
        else
            β₁,β₂,β₃ = Δ₁,Δ₂,ψ;
        end;

        # 2 compute the scattering coefficients
        so = cosd(sza) * cosd(vza);
        T₁ = 2 * Cs * Co + Ss * So * cos(ψ);
        T₂ = sin(β₂) * (2 * Ds * Do + Ss * So * cos(β₁) * cos(β₃));
        F₁ = ((FT(π) - β₂) * T₁ + T₂) / (2 * so * FT(π));
        F₂ = (-β₂ * T₁ + T₂) / (2 * so * FT(π));

        # 3 compute the area scattering coefficient fractions (sb for backward and sf for forward)
        sensa.sb_incl[i] = (F₂ >= 0 ? F₁ : abs(F₂));
        sensa.sf_incl[i] = (F₂ >= 0 ? F₂ : abs(F₁));
    end;
    sensa.ko_leaf = cansa.p_incl_leaf' * sensa.ko_incl * sensa.ci_sensor;
    sensa.ko_stem = cansa.p_incl_stem' * sensa.ko_incl * sensa.ci_sensor;

    # compute the scattering weights for diffuse/direct -> sensor for backward and forward scattering
    sensa.w_dob_leaf = 0;
    sensa.w_dof_leaf = 0;
    sensa.w_dob_stem = 0;
    sensa.w_dof_stem = 0;
    for i in eachindex(Θ_INCL)
        f_ada = f_adaxial(senst.vza, Θ_INCL[i]);
        f_aba = 1 - f_ada;
        f_inc = Θ_INCL[i] / 180;
        sensa.w_dob_leaf += (f_ada * (1 - f_inc) + f_aba * f_inc) * cansa.p_incl_leaf[i] * sensa.ko_incl[i];
        sensa.w_dof_leaf += (f_ada * f_inc + f_aba * (1 - f_inc)) * cansa.p_incl_leaf[i] * sensa.ko_incl[i];
        sensa.w_dob_stem += (f_ada * (1 - f_inc) + f_aba * f_inc) * cansa.p_incl_stem[i] * sensa.ko_incl[i];
        sensa.w_dof_stem += (f_ada * f_inc + f_aba * (1 - f_inc)) * cansa.p_incl_stem[i] * sensa.ko_incl[i];
    end;
    sensa.w_sob_leaf = cansa.p_incl_leaf' * sensa.sb_incl;
    sensa.w_sof_leaf = cansa.p_incl_leaf' * sensa.sf_incl;
    sensa.w_sob_stem = cansa.p_incl_stem' * sensa.sb_incl;
    sensa.w_sof_stem = cansa.p_incl_stem' * sensa.sf_incl;

    # compute the fo and fo_abs matrices
    for i in eachindex(Θ_AZI)
        cos_azi_raa = cosd(Θ_AZI[i] .- (senst.vaa - sunst.saa));
        @. sensa.fo[:,i] = sensa.Co_incl + sensa.So_incl * cos_azi_raa;
    end;
    @. sensa.fo /= cosd(senst.vza);
    @. sensa.fo_abs = abs(sensa.fo);
    for i in eachindex(Θ_INCL)
        @. sensa.fo_cos²_incl[i,:] = (@view sensa.fo[i,:]) * cosd(Θ_INCL[i]) ^ 2; # TODO: is this related to the bf calculation in SCOPE?
    end;
    @. sensa.fo_fs = sunsa.fs * sensa.fo;
    @. sensa.fo_fs_abs = abs(sensa.fo_fs);

    # compute fractions of leaves/soil that can be viewed from the sensor direction
    #     it is different from the SCOPE model that we compute the po directly for canopy layers rather than the boundaries (last one is still soil though)
    kocipai = sensa.ko_leaf * canst.lai + sensa.ko_stem * canst.sai;
    for i in eachindex(canst.δlai)
        kociipai = sensa.ko_leaf * canst.δlai[i] + sensa.ko_stem * canst.δsai[i];
        sensa.p_sensor[i] = sensa.ci_sensor / kociipai * (exp(kocipai * cansa.x_bnds[i]) - exp(kocipai * cansa.x_bnds[i+1]));
    end;
    sensa.p_sensor_soil = exp(-kocipai);

    # TODO: the pso function could lead to pso > ps or pso > po, redo the calculation without using the min function
    # compute the fraction of sunlit leaves that can be viewed from the sensor direction (for hot spot)
    # equations from Appendix C of the mSCOPE paper (Yang et al., 2018)
    pai = canst.lai + canst.sai;
    ag = sqrt( tand(sunst.sza) ^ 2 + tand(senst.vza) ^ 2 - 2 * tand(sunst.sza) * tand(senst.vza) * cosd(senst.vaa - sunst.saa) );
    Σk = (sunsa.ks_leaf * canst.lai + sunsa.ks_stem * canst.sai + sensa.ko_leaf * canst.lai + sensa.ko_stem * canst.sai);
    Πk = sqrt((sunsa.ks_leaf * canst.lai + sunsa.ks_stem * canst.sai) * (sensa.ko_leaf * canst.lai + sensa.ko_stem * canst.sai));
    sl = lw2ch * 2 * pai / Σk;
    pso(x) = ag == 0 ? sensa.ci_sensor * exp(Σk * x - Πk * x) : sensa.ci_sensor * exp(Σk * x + Πk * sl / ag * (1 - exp(ag / sl * x)));

    for i in eachindex(canst.δlai)
        sensa.p_sun_sensor[i] = quadgk(pso, cansa.x_bnds[i+1], cansa.x_bnds[i]; rtol = 1e-4)[1] / (cansa.x_bnds[i] - cansa.x_bnds[i+1]);
        sensa.p_sun_sensor[i] = min(sensa.p_sun_sensor[i], sensa.p_sensor[i], sunsa.p_sunlit[i]);
    end;

    return nothing
);

sensor_geometry_aux!(
            config::SPACConfig{FT},
            ::CanopyRTSCOPE,
            canst::CanopyStructureTrait{FT},
            cansa::CanopyStructureAuxil{FT},
            sunst::SunGeometryState{FT},
            sunsa::SunGeometryAuxil{FT},
            senst::SensorGeometryState{FT},
            sensa::SensorGeometryAuxil{FT},
            lw2ch::FT) where {FT} = (
    # if none of REF or SIF is enabled, or sza > 89, or LAI+SAI <= 0, do nothing
    if (!config.FEATURES.ENABLE_REF && !config.FEATURES.ENABLE_SIF) || sunst.sza > 89 || (canst.lai <= 0 && canst.sai <= 0)
        return nothing
    end;

    (; Θ_AZI, Θ_INCL) = config.DIMENSIONS;

    # compute clumping index from sensor zenith angle
    sensa.ci_sensor = canst.ci.ci_0 * (1 - canst.ci.ci_1 * cosd(senst.vza));

    # extinction coefficients for the solar radiation
    vza = senst.vza;
    sza = sunst.sza;
    raa = senst.vaa - sunst.saa;
    for i in eachindex(Θ_INCL)
        # observer angle
        Co = cosd(Θ_INCL[i]) * cosd(vza);
        So = sind(Θ_INCL[i]) * sind(vza);
        cosβo = abs(So) <= 1e-6 ? FT(1) : -Co/So;
        if abs(cosβo) < 1
            βo = acos(cosβo);
            Do = So;
        elseif vza < 90
            βo = FT(π);
            Do = Co;
        else
            βo = 0;
            Do = -Co;
        end;
        sensa.Co_incl[i] = Co;
        sensa.So_incl[i] = So;
        sensa.βo_incl[i] = βo;
        sensa.ko_incl[i] = 2 / FT(π) / cosd(vza) * (Co * (βo - FT(π)/2) + So * sin(βo));

        # compute the scattering coefficients
        Cs = sunsa.Cs_incl[i];
        Ss = sunsa.Ss_incl[i];
        βs = sunsa.βs_incl[i];
        Ds = (abs(cos(βs)) < 1 ? Ss : Cs);

        # 1 compute the Δ and β angles
        Δ₁ = abs(βs - βo);
        Δ₂ = FT(π) - abs(βs + βo - FT(π));

        ψ = deg2rad( abs(raa - 360*round(raa/360)) );
        if ψ <= Δ₁
            β₁,β₂,β₃ = ψ,Δ₁,Δ₂;
        elseif Δ₁ < ψ < Δ₂
            β₁,β₂,β₃ = Δ₁,ψ,Δ₂;
        else
            β₁,β₂,β₃ = Δ₁,Δ₂,ψ;
        end;

        # 2 compute the scattering coefficients
        so = cosd(sza) * cosd(vza);
        T₁ = 2 * Cs * Co + Ss * So * cos(ψ);
        T₂ = sin(β₂) * (2 * Ds * Do + Ss * So * cos(β₁) * cos(β₃));
        F₁ = ((FT(π) - β₂) * T₁ + T₂) / (2 * so * FT(π));
        F₂ = (-β₂ * T₁ + T₂) / (2 * so * FT(π));

        # 3 compute the area scattering coefficient fractions (sb for backward and sf for forward)
        sensa.sb_incl[i] = (F₂ >= 0 ? F₁ : abs(F₂));
        sensa.sf_incl[i] = (F₂ >= 0 ? F₂ : abs(F₁));
    end;
    sensa.ko_leaf = cansa.p_incl_leaf' * sensa.ko_incl * sensa.ci_sensor;
    sensa.ko_stem = cansa.p_incl_stem' * sensa.ko_incl * sensa.ci_sensor;

    # compute the scattering weights for diffuse/direct -> sensor for backward and forward scattering
    sensa.w_dob_leaf = (sensa.ko_leaf + cansa.bf_leaf) / 2;
    sensa.w_dof_leaf = (sensa.ko_leaf - cansa.bf_leaf) / 2;
    sensa.w_dob_stem = (sensa.ko_stem + cansa.bf_stem) / 2;
    sensa.w_dof_stem = (sensa.ko_stem - cansa.bf_stem) / 2;
    sensa.w_sob_leaf = cansa.p_incl_leaf' * sensa.sb_incl;
    sensa.w_sof_leaf = cansa.p_incl_leaf' * sensa.sf_incl;
    sensa.w_sob_stem = cansa.p_incl_stem' * sensa.sb_incl;
    sensa.w_sof_stem = cansa.p_incl_stem' * sensa.sf_incl;

    # compute the fo and fo_abs matrices
    for i in eachindex(Θ_AZI)
        cos_azi_raa = cosd(Θ_AZI[i] .- (senst.vaa - sunst.saa));
        @. sensa.fo[:,i] = sensa.Co_incl + sensa.So_incl * cos_azi_raa;
    end;
    @. sensa.fo /= cosd(senst.vza);
    @. sensa.fo_abs = abs(sensa.fo);
    for i in eachindex(Θ_INCL)
        @. sensa.fo_cos²_incl[i,:] = (@view sensa.fo[i,:]) * cosd(Θ_INCL[i]) ^ 2; # TODO: is this related to the bf calculation in SCOPE?
    end;
    @. sensa.fo_fs = sunsa.fs * sensa.fo;
    @. sensa.fo_fs_abs = abs(sensa.fo_fs);

    # compute fractions of leaves/soil that can be viewed from the sensor direction
    #     it is different from the SCOPE model that we compute the po directly for canopy layers rather than the boundaries (last one is still soil though)
    kocipai = sensa.ko_leaf * canst.lai + sensa.ko_stem * canst.sai;
    for i in eachindex(canst.δlai)
        kociipai = sensa.ko_leaf * canst.δlai[i] + sensa.ko_stem * canst.δsai[i];
        sensa.p_sensor[i] = sensa.ci_sensor / kociipai * (exp(kocipai * cansa.x_bnds[i]) - exp(kocipai * cansa.x_bnds[i+1]));
    end;
    sensa.p_sensor_soil = exp(-kocipai);

    # TODO: the pso function could lead to pso > ps or pso > po, redo the calculation without using the min function
    # compute the fraction of sunlit leaves that can be viewed from the sensor direction (for hot spot)
    # equations from Appendix C of the mSCOPE paper (Yang et al., 2018)
    pai = canst.lai + canst.sai;
    ag = sqrt( tand(sunst.sza) ^ 2 + tand(senst.vza) ^ 2 - 2 * tand(sunst.sza) * tand(senst.vza) * cosd(senst.vaa - sunst.saa) );
    Σk = (sunsa.ks_leaf * canst.lai + sunsa.ks_stem * canst.sai + sensa.ko_leaf * canst.lai + sensa.ko_stem * canst.sai);
    Πk = sqrt((sunsa.ks_leaf * canst.lai + sunsa.ks_stem * canst.sai) * (sensa.ko_leaf * canst.lai + sensa.ko_stem * canst.sai));
    sl = lw2ch * 2 * pai / Σk;
    pso(x) = ag == 0 ? sensa.ci_sensor * exp(Σk * x - Πk * x) : sensa.ci_sensor * exp(Σk * x + Πk * sl / ag * (1 - exp(ag / sl * x)));

    for i in eachindex(canst.δlai)
        sensa.p_sun_sensor[i] = quadgk(pso, cansa.x_bnds[i+1], cansa.x_bnds[i]; rtol = 1e-4)[1] / (cansa.x_bnds[i] - cansa.x_bnds[i+1]);
        sensa.p_sun_sensor[i] = min(sensa.p_sun_sensor[i], sensa.p_sensor[i], sunsa.p_sunlit[i]);
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sensor_geometry! (run per viewing zenith angle)
#     2023-Oct-11: compute canopy layer scattering
#     2023-Oct-13: improve p_sun_sensor calculation accuracy
#     2023-Oct-14: if none of REF or SIF is enabled, skip the sensor geometry calculation
#     2023-Oct-18: account for SAI in the sensor geometry calculation
#     2024-Feb-22: add solar zenith angle control
#     2024-Feb-25: move the trait- and state-dependent calculations to the sensor_geometry_aux! function
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Oct-16: add option to compute effective leaf spectra based on CI
#     2024-Nov-08: when using EFFECTIVE_LEAF_SPECTRA make sure LAI > 0
#
#######################################################################################################################################################################################################
"""

    sensor_geometry!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sensor_geometry!(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;
    sun_geo = spac.canopy.sun_geometry;

    if (!config.FEATURES.ENABLE_REF && !config.FEATURES.ENABLE_SIF) || sun_geo.state.sza > 89 || (can_str.trait.lai <= 0 && can_str.trait.sai <= 0)
        return nothing
    end;

    # run the sensor geometry simulations only if any of canopy reflectance feature or fluorescence feature is enabled and if LAI+SAI > 0
    (; SPECTRA) = config.CONSTANTS;
    (; EFFECTIVE_LEAF_SPECTRA) = config.FEATURES;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    n_layer = length(leaves);

    # compute the effective leaf reflectance and transmittance spectra using the PROSPECT model scheme
    mask_effective = EFFECTIVE_LEAF_SPECTRA && can_str.trait.lai > 0;
    if mask_effective
        ρ_2 = spac.cache.cache_wl_1;
        τ_2 = spac.cache.cache_wl_2;
        n_eff = 1 / sen_geo.auxil.ci_sensor;
        for irt in 1:n_layer
            ilf = n_layer + 1 - irt;
            leaf = leaves[ilf];
            ρ_1 = leaf.bio.auxil.ρ_leaf;
            τ_1 = leaf.bio.auxil.τ_leaf;
            ρ_2 .= layer_2_ρ.(ρ_1, τ_1, n_eff - 1);
            τ_2 .= layer_2_τ.(ρ_1, τ_1, n_eff - 1);
            sen_geo.auxil.ρ_leaf_eff[:,irt] .= leaf_ρ.(ρ_1, τ_1, ρ_1, τ_1, ρ_2);
            sen_geo.auxil.τ_leaf_eff[:,irt] .= leaf_τ.(τ_1, ρ_1, ρ_2, τ_2);
        end;
    end;

    # compute the scattering coefficients per leaf area
    ρ_leaf_dif = spac.cache.cache_wl_1;
    τ_leaf_dif = spac.cache.cache_wl_2;
    ρ_leaf_dir = spac.cache.cache_wl_3;
    τ_leaf_dir = spac.cache.cache_wl_4;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        mask_effective ? ρ_leaf_dif .= view(sen_geo.auxil.ρ_leaf_eff,:,irt) : ρ_leaf_dif .= leaf.bio.auxil.ρ_leaf;
        mask_effective ? τ_leaf_dif .= view(sen_geo.auxil.τ_leaf_eff,:,irt) : τ_leaf_dif .= leaf.bio.auxil.τ_leaf;
        mask_effective ? ρ_leaf_dir .= view(sun_geo.auxil.ρ_leaf_eff,:,irt) : ρ_leaf_dir .= leaf.bio.auxil.ρ_leaf;
        mask_effective ? τ_leaf_dir .= view(sun_geo.auxil.τ_leaf_eff,:,irt) : τ_leaf_dir .= leaf.bio.auxil.τ_leaf;
        @. sen_geo.auxil.dob_leaf[:,irt] = sen_geo.auxil.w_dob_leaf * ρ_leaf_dif + sen_geo.auxil.w_dof_leaf * τ_leaf_dif;
        @. sen_geo.auxil.dof_leaf[:,irt] = sen_geo.auxil.w_dof_leaf * ρ_leaf_dif + sen_geo.auxil.w_dob_leaf * τ_leaf_dif;
        @. sen_geo.auxil.so_leaf[:,irt]  = sen_geo.auxil.w_sob_leaf * ρ_leaf_dir + sen_geo.auxil.w_sof_leaf * τ_leaf_dir;
        @. sen_geo.auxil.dob_stem[:,irt] = sen_geo.auxil.w_dob_stem * SPECTRA.ρ_STEM;
        @. sen_geo.auxil.dof_stem[:,irt] = sen_geo.auxil.w_dof_stem * SPECTRA.ρ_STEM;
        @. sen_geo.auxil.so_stem[:,irt]  = sen_geo.auxil.w_sob_stem * SPECTRA.ρ_STEM;
    end;

    return nothing
end;

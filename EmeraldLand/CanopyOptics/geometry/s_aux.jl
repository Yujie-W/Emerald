#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-25: add s_aux! method for SunGeometryState-dependent variables
#     2024-Feb-25: add s_aux! method for SensorGeometryState-dependent variables
#     2024-Feb-25: add s_aux! method for the combined MultiLayerCanopy
#
#######################################################################################################################################################################################################
s_aux!(config::SPACConfiguration{FT}, trait::CanopyStructureTrait{FT}, t_aux::CanopyStructureTDAuxil{FT}, sunst::SunGeometryState{FT}, sunsa::SunGeometrySDAuxil{FT}) where {FT} = (
    # if sza > 89 or both LAI and SAI are zero, do nothing
    if sunst.sza > 89 || (trait.lai <= 0 && trait.sai <= 0)
        return nothing
    end;

    (; Θ_AZI, Θ_INCL) = config;

    # extinction coefficients for the solar radiation
    for i in eachindex(Θ_INCL)
        Cs = cosd(Θ_INCL[i]) * cosd(sunst.sza);
        Ss = sind(Θ_INCL[i]) * sind(sunst.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        sunsa.Cs_incl[i] = Cs;
        sunsa.Ss_incl[i] = Ss;
        sunsa.βs_incl[i] = βs;
        sunsa.ks_incl[i] = 2 / FT(π) / cosd(sunst.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;
    sunsa.ks = t_aux.p_incl' * sunsa.ks_incl;

    # compute the scattering weights for diffuse/direct -> diffuse for backward and forward scattering
    sunsa.sdb = (sunsa.ks + t_aux.bf) / 2;
    sunsa.sdf = (sunsa.ks - t_aux.bf) / 2;

    # compute the sunlit leaf fraction
    # sunsa.ps = exp.(sunsa.ks .* trait.ci * trait.lai .* t_aux.x_bnds);
    kscipai = sunsa.ks * trait.ci * (trait.lai + trait.sai);
    for i in eachindex(trait.δlai)
        ksciipai = sunsa.ks * trait.ci * (trait.δlai[i] + trait.δsai[i]);
        sunsa.p_sunlit[i] = 1 / ksciipai * (exp(kscipai * t_aux.x_bnds[i]) - exp(kscipai * t_aux.x_bnds[i+1]));
    end;

    # compute the fs and fs_abs matrices
    for i in eachindex(Θ_AZI)
        view(sunsa.fs,:,i) .= sunsa.Cs_incl .+ sunsa.Ss_incl .* cosd(Θ_AZI[i]);
    end;
    sunsa.fs ./= cosd(sunst.sza);
    sunsa.fs_abs .= abs.(sunsa.fs);
    mul!(sunsa.fs_abs_mean, sunsa.fs_abs', t_aux.p_incl);
    for i in eachindex(Θ_INCL)
        view(sunsa.fs_cos²_incl,i,:) .= view(sunsa.fs,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;

    return nothing
);

s_aux!(config::SPACConfiguration{FT},
       trait::CanopyStructureTrait{FT},
       t_aux::CanopyStructureTDAuxil{FT},
       sunst::SunGeometryState{FT},
       sunsa::SunGeometrySDAuxil{FT},
       senst::SensorGeometryState{FT},
       sensa::SensorGeometrySDAuxil{FT}) where {FT} = (
    # if none of REF or SIF is enabled, or sza > 89, or LAI+SAI <= 0, do nothing
    if (!config.ENABLE_REF && !config.ENABLE_SIF) || sunst.sza > 89 || (trait.lai <= 0 && trait.sai <= 0)
        return nothing
    end;

    (; Θ_AZI, Θ_INCL) = config;

    # extinction coefficients for the solar radiation
    vza = senst.vza;
    sza = sunst.sza;
    raa = senst.vaa - sunst.saa;
    for i in eachindex(Θ_INCL)
        Co = cosd(Θ_INCL[i]) * cosd(vza);
        So = sind(Θ_INCL[i]) * sind(vza);
        βo = (Co >= So ? FT(π) : acos(-Co/So));
        sensa.Co_incl[i] = Co;
        sensa.So_incl[i] = So;
        sensa.βo_incl[i] = βo;
        sensa.ko_incl[i] = 2 / FT(π) / cosd(FT(π)) * (Co * (βo - FT(π)/2) + So * sin(βo));

        # compute the scattering coefficients
        Cs = sunsa.Cs_incl[i];
        Ss = sunsa.Ss_incl[i];
        βs = sunsa.βs_incl[i];

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
        Ds = (βs < pi ? Ss : Cs);
        Do = (0 < βo < pi ? So : -Co/cos(βo));
        so = cosd(sza) * cosd(vza);
        T₁ = 2 * Cs * Co + Ss * So * cosd(ψ);
        T₂ = sin(β₂) * (2 * Ds * Do + Ss * So * cos(β₁) * cos(β₃));
        F₁ = ((FT(π) - β₂) * T₁ + T₂) / abs(so);
        F₂ = (-β₂ * T₁ + T₂) / abs(so);

        # 3 compute the area scattering coefficient fractions (sb for backward and sf for forward)
        sensa.sb_incl[i] = (F₂ >= 0 ? F₁ : abs(F₂)) / (2 * FT(π));
        sensa.sf_incl[i] = (F₂ >= 0 ? F₂ : abs(F₁)) / (2 * FT(π));
    end;
    sensa.ko = t_aux.p_incl' * sensa.ko_incl;

    # compute the scattering weights for diffuse/direct -> sensor for backward and forward scattering
    sensa.dob = (sensa.ko + t_aux.bf) / 2;
    sensa.dof = (sensa.ko - t_aux.bf) / 2;
    sensa.sob = t_aux.p_incl' * sensa.sb_incl;
    sensa.sof = t_aux.p_incl' * sensa.sf_incl;

    # compute the fo and fo_abs matrices
    for i in eachindex(Θ_AZI)
        cos_azi_raa = cosd(Θ_AZI[i] .- (senst.vaa - sunst.saa));
        view(sensa.fo,:,i) .= sensa.Co_incl .+ sensa.So_incl .* cos_azi_raa;
    end;
    sensa.fo ./= cosd(senst.vza);
    sensa.fo_abs .= abs.(sensa.fo);
    for i in eachindex(Θ_INCL)
        view(sensa.fo_cos²_incl,i,:) .= view(sensa.fo,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;
    sensa.fo_fs .= sunsa.fs .* sensa.fo;
    sensa.fo_fs_abs .= abs.(sensa.fo_fs);

    # compute fractions of leaves/soil that can be viewed from the sensor direction
    #     it is different from the SCOPE model that we compute the po directly for canopy layers rather than the boundaries (last one is still soil though)
    kocipai = sensa.ko * trait.ci * (trait.lai + trait.sai);
    for i in eachindex(trait.δlai)
        kociipai = sensa.ko * trait.ci * (trait.δlai[i] + trait.δsai[i]);
        sensa.p_sensor[i] = 1 / kociipai * (exp(kocipai * t_aux.x_bnds[i]) - exp(kocipai * t_aux.x_bnds[i+1]));
    end;
    sensa.p_sensor_soil = exp(-kocipai);

    # compute the fraction of sunlit leaves that can be viewed from the sensor direction (for hot spot)
    dso = sqrt( tand(sunst.sza) ^ 2 + tand(senst.vza) ^ 2 - 2 * tand(sunst.sza) * tand(senst.vza) * cosd(senst.vaa - sunst.saa) );
    Σk = sensa.ko + sunsa.ks;
    Πk = sensa.ko * sunsa.ks;
    cl = trait.ci * (trait.lai + trait.sai);
    α  = dso / trait.hot_spot * 2 / Σk;
    pso(x) = dso == 0 ? exp( (Σk - sqrt(Πk)) * cl * x ) : exp( Σk * cl * x + sqrt(Πk) * cl / α * (1 - exp(α * x)) );

    for i in eachindex(trait.δlai)
        sensa.p_sun_sensor[i] = quadgk(pso, t_aux.x_bnds[i+1], t_aux.x_bnds[i]; rtol = 1e-2)[1] / (t_aux.x_bnds[i] - t_aux.x_bnds[i+1]);
    end;

    return nothing
);

s_aux!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT} = (
    s_aux!(config, can.structure.trait, can.structure.t_aux, can.sun_geometry.state, can.sun_geometry.s_aux);
    s_aux!(config, can.structure.trait, can.structure.t_aux, can.sun_geometry.state, can.sun_geometry.s_aux, can.sensor_geometry.state, can.sensor_geometry.s_aux);

    return nothing
);

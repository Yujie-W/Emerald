# This file contains functions to compute the sensor geometry of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-22: add function sensor_geometry_aux! to update the state-dependent auxiliary variables for sensor geometry (to call in step_remote_sensing!)
#
#######################################################################################################################################################################################################
"""

    sensor_geometry_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sensor_geometry_aux! end;

sensor_geometry_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = sensor_geometry_aux!(config, spac.canopy);

sensor_geometry_aux!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT} = (
    sensor_geometry_aux!(config, can.structure.trait, can.structure.t_aux, can.sun_geometry.state, can.sun_geometry.s_aux, can.sensor_geometry.state, can.sensor_geometry.s_aux);

    return nothing
);

sensor_geometry_aux!(
            config::SPACConfiguration{FT},
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

    (; HOT_SPOT, Θ_AZI, Θ_INCL) = config;

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
    α  = dso / HOT_SPOT * 2 / Σk;
    pso(x) = dso == 0 ? exp( (Σk - sqrt(Πk)) * cl * x ) : exp( Σk * cl * x + sqrt(Πk) * cl / α * (1 - exp(α * x)) );

    for i in eachindex(trait.δlai)
        sensa.p_sun_sensor[i] = quadgk(pso, t_aux.x_bnds[i+1], t_aux.x_bnds[i]; rtol = 1e-2)[1] / (t_aux.x_bnds[i] - t_aux.x_bnds[i+1]);
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
#
#######################################################################################################################################################################################################
"""

    sensor_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sensor_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;
    sun_geo = spac.canopy.sun_geometry;

    if (!config.ENABLE_REF && !config.ENABLE_SIF) || sun_geo.state.sza > 89 || (can_str.trait.lai <= 0 && can_str.trait.sai <= 0)
        return nothing
    end;

    # run the sensor geometry simulations only if any of canopy reflectance feature or fluorescence feature is enabled and if LAI+SAI > 0
    (; SPECTRA) = config;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    n_layer = length(leaves);

    # compute the scattering coefficients per leaf area
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sen_geo.auxil.dob_leaf[:,irt] .= sen_geo.s_aux.dob * leaf.bio.auxil.ρ_leaf .+ sen_geo.s_aux.dof * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.dof_leaf[:,irt] .= sen_geo.s_aux.dof * leaf.bio.auxil.ρ_leaf .+ sen_geo.s_aux.dob * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.so_leaf[:,irt]  .= sen_geo.s_aux.sob * leaf.bio.auxil.ρ_leaf .+ sen_geo.s_aux.sof * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.dob_stem[:,irt] .= sen_geo.s_aux.dob * SPECTRA.ρ_STEM;
        sen_geo.auxil.dof_stem[:,irt] .= sen_geo.s_aux.dof * SPECTRA.ρ_STEM;
        sen_geo.auxil.so_stem[:,irt]  .= sen_geo.s_aux.sob * SPECTRA.ρ_STEM;
    end;

    return nothing
end;

# This file contains functions to compute the sensor geometry of the canopy

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
    (; SPECTRA, Θ_AZI, Θ_INCL) = config;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    n_layer = length(leaves);

    # extinction coefficients for the solar radiation
    vza = sen_geo.state.vza;
    sza = sun_geo.state.sza;
    raa = sen_geo.state.vaa - sun_geo.state.saa;
    for i in eachindex(Θ_INCL)
        Co = cosd(Θ_INCL[i]) * cosd(vza);
        So = sind(Θ_INCL[i]) * sind(vza);
        βo = (Co >= So ? FT(π) : acos(-Co/So));
        sen_geo.auxil.Co_incl[i] = Co;
        sen_geo.auxil.So_incl[i] = So;
        sen_geo.auxil.βo_incl[i] = βo;
        sen_geo.auxil.ko_incl[i] = 2 / FT(π) / cosd(FT(π)) * (Co * (βo - FT(π)/2) + So * sin(βo));

        # compute the scattering coefficients
        Cs = sun_geo.s_aux.Cs_incl[i];
        Ss = sun_geo.s_aux.Ss_incl[i];
        βs = sun_geo.s_aux.βs_incl[i];

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
        sen_geo.auxil.sb_incl[i] = (F₂ >= 0 ? F₁ : abs(F₂)) / (2 * FT(π));
        sen_geo.auxil.sf_incl[i] = (F₂ >= 0 ? F₂ : abs(F₁)) / (2 * FT(π));
    end;
    sen_geo.auxil.ko = can_str.t_aux.p_incl' * sen_geo.auxil.ko_incl;

    # compute the scattering weights for diffuse/direct -> sensor for backward and forward scattering
    sen_geo.auxil.dob = (sen_geo.auxil.ko + can_str.t_aux.bf) / 2;
    sen_geo.auxil.dof = (sen_geo.auxil.ko - can_str.t_aux.bf) / 2;
    sen_geo.auxil.sob = can_str.t_aux.p_incl' * sen_geo.auxil.sb_incl;
    sen_geo.auxil.sof = can_str.t_aux.p_incl' * sen_geo.auxil.sf_incl;

    # compute the fo and fo_abs matrices
    for i in eachindex(Θ_AZI)
        cos_azi_raa = cosd(Θ_AZI[i] .- (sen_geo.state.vaa - sun_geo.state.saa));
        view(sen_geo.auxil.fo,:,i) .= sen_geo.auxil.Co_incl .+ sen_geo.auxil.So_incl .* cos_azi_raa;
    end;
    sen_geo.auxil.fo ./= cosd(sen_geo.state.vza);
    sen_geo.auxil.fo_abs .= abs.(sen_geo.auxil.fo);
    for i in eachindex(Θ_INCL)
        view(sen_geo.auxil.fo_cos²_incl,i,:) .= view(sen_geo.auxil.fo,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;
    sen_geo.auxil.fo_fs .= sun_geo.s_aux.fs .* sen_geo.auxil.fo;
    sen_geo.auxil.fo_fs_abs .= abs.(sen_geo.auxil.fo_fs);

    # compute fractions of leaves/soil that can be viewed from the sensor direction
    #     it is different from the SCOPE model that we compute the po directly for canopy layers rather than the boundaries (last one is still soil though)
    kocipai = sen_geo.auxil.ko * can_str.trait.ci * (can_str.trait.lai + can_str.trait.sai);
    for i in 1:n_layer
        kociipai = sen_geo.auxil.ko * can_str.trait.ci * (can_str.trait.δlai[i] + can_str.trait.δsai[i]);
        sen_geo.auxil.p_sensor[i] = 1 / kociipai * (exp(kocipai * can_str.t_aux.x_bnds[i]) - exp(kocipai * can_str.t_aux.x_bnds[i+1]));
    end;
    sen_geo.auxil.p_sensor_soil = exp(-kocipai);

    # compute the fraction of sunlit leaves that can be viewed from the sensor direction (for hot spot)
    dso = sqrt( tand(sun_geo.state.sza) ^ 2 + tand(sen_geo.state.vza) ^ 2 - 2 * tand(sun_geo.state.sza) * tand(sen_geo.state.vza) * cosd(sen_geo.state.vaa - sun_geo.state.saa) );
    Σk = sen_geo.auxil.ko + sun_geo.s_aux.ks;
    Πk = sen_geo.auxil.ko * sun_geo.s_aux.ks;
    cl = can_str.trait.ci * (can_str.trait.lai + can_str.trait.sai);
    α  = dso / can_str.trait.hot_spot * 2 / Σk;
    pso(x) = dso == 0 ? exp( (Σk - sqrt(Πk)) * cl * x ) : exp( Σk * cl * x + sqrt(Πk) * cl / α * (1 - exp(α * x)) );

    for i in 1:n_layer
        sen_geo.auxil.p_sun_sensor[i] = quadgk(pso, can_str.t_aux.x_bnds[i+1], can_str.t_aux.x_bnds[i]; rtol = 1e-2)[1] / (can_str.t_aux.x_bnds[i] - can_str.t_aux.x_bnds[i+1]);
    end;

    # compute the scattering coefficients per leaf area
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sen_geo.auxil.dob_leaf[:,irt] .= sen_geo.auxil.dob * leaf.bio.auxil.ρ_leaf .+ sen_geo.auxil.dof * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.dof_leaf[:,irt] .= sen_geo.auxil.dof * leaf.bio.auxil.ρ_leaf .+ sen_geo.auxil.dob * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.so_leaf[:,irt]  .= sen_geo.auxil.sob * leaf.bio.auxil.ρ_leaf .+ sen_geo.auxil.sof * leaf.bio.auxil.τ_leaf;
        sen_geo.auxil.dob_stem[:,irt] .= sen_geo.auxil.dob * SPECTRA.ρ_STEM;
        sen_geo.auxil.dof_stem[:,irt] .= sen_geo.auxil.dof * SPECTRA.ρ_STEM;
        sen_geo.auxil.so_stem[:,irt]  .= sen_geo.auxil.sob * SPECTRA.ρ_STEM;
    end;

    return nothing
end;

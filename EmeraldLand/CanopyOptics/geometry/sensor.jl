# This file contains functions to compute the sensor geometry of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sensor_geometry! (run per viewing zenith angle)
#     2023-Oct-11: compute canopy layer scattering
#     2023-Oct-13: improve p_sun_sensor calculation accuracy
#
#######################################################################################################################################################################################################
"""

    sensor_geometry!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sensor_geometry!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; DIM_LAYER, Θ_AZI, Θ_INCL) = config;
    (; CANOPY, LEAVES) = spac;

    # extinction coefficients for the solar radiation
    vza = CANOPY.sensor_geometry.state.vza;
    sza = CANOPY.sun_geometry.state.sza;
    raa = CANOPY.sensor_geometry.state.vaa - CANOPY.sun_geometry.state.saa;
    for i in eachindex(Θ_INCL)
        Co = cosd(Θ_INCL[i]) * cosd(vza);
        So = sind(Θ_INCL[i]) * sind(vza);
        βo = (Co >= So ? FT(π) : acos(-Co/So));
        CANOPY.sensor_geometry.auxil.Co_incl[i] = Co;
        CANOPY.sensor_geometry.auxil.So_incl[i] = So;
        CANOPY.sensor_geometry.auxil.βo_incl[i] = βo;
        CANOPY.sensor_geometry.auxil.ko_incl[i] = 2 / FT(π) / cosd(FT(π)) * (Co * (βo - FT(π)/2) + So * sin(βo));

        # compute the scattering coefficients
        Cs = CANOPY.sun_geometry.auxil.Cs_incl[i];
        Ss = CANOPY.sun_geometry.auxil.Ss_incl[i];
        βs = CANOPY.sun_geometry.auxil.βs_incl[i];

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
        CANOPY.sensor_geometry.auxil.sb_incl[i] = (F₂ >= 0 ? F₁ : abs(F₂)) / (2 * FT(π));
        CANOPY.sensor_geometry.auxil.sf_incl[i] = (F₂ >= 0 ? F₂ : abs(F₁)) / (2 * FT(π));
    end;
    CANOPY.sensor_geometry.auxil.ko = CANOPY.structure.state.p_incl' * CANOPY.sensor_geometry.auxil.ko_incl;

    # compute the scattering weights for diffuse/direct -> sensor for backward and forward scattering
    CANOPY.sensor_geometry.auxil.dob = (CANOPY.sensor_geometry.auxil.ko + CANOPY.structure.auxil.bf) / 2;
    CANOPY.sensor_geometry.auxil.dof = (CANOPY.sensor_geometry.auxil.ko - CANOPY.structure.auxil.bf) / 2;
    CANOPY.sensor_geometry.auxil.sob = CANOPY.structure.state.p_incl' * CANOPY.sensor_geometry.auxil.sb_incl;
    CANOPY.sensor_geometry.auxil.sof = CANOPY.structure.state.p_incl' * CANOPY.sensor_geometry.auxil.sf_incl;

    # compute the fo and fo_abs matrices
    for i in eachindex(Θ_AZI)
        cos_azi_raa = cosd(Θ_AZI[i] .- (CANOPY.sensor_geometry.state.vaa - CANOPY.sun_geometry.state.saa));
        view(CANOPY.sensor_geometry.auxil.fo,:,i) .= CANOPY.sensor_geometry.auxil.Co_incl .+ CANOPY.sensor_geometry.auxil.So_incl .* cos_azi_raa;
    end;
    CANOPY.sensor_geometry.auxil.fo ./= cosd(CANOPY.sensor_geometry.state.vza);
    CANOPY.sensor_geometry.auxil.fo_abs .= abs.(CANOPY.sensor_geometry.auxil.fo);
    for i in eachindex(Θ_INCL)
        view(CANOPY.sensor_geometry.auxil.fo_cos²_incl,i,:) .= view(CANOPY.sensor_geometry.auxil.fo,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;
    CANOPY.sensor_geometry.auxil.fo_fs .= CANOPY.sun_geometry.auxil.fs .* CANOPY.sensor_geometry.auxil.fo;
    CANOPY.sensor_geometry.auxil.fo_fs_abs .= abs.(CANOPY.sensor_geometry.auxil.fo_fs);

    # compute fractions of leaves/soil that can be viewed from the sensor direction
    #     it is different from the SCOPE model that we compute the po directly for canopy layers rather than the boundaries (last one is still soil though)
    kocilai = CANOPY.sensor_geometry.auxil.ko * CANOPY.structure.auxil.ci * CANOPY.structure.state.lai;
    for i in 1:DIM_LAYER
        koilai = CANOPY.sensor_geometry.auxil.ko * CANOPY.structure.auxil.ci * CANOPY.structure.state.δlai[i];
        CANOPY.sensor_geometry.auxil.p_sensor[i] = 1 / koilai * (exp(kocilai * CANOPY.structure.auxil.x_bnds[i]) - exp(kocilai * CANOPY.structure.auxil.x_bnds[i+1]));
    end;
    CANOPY.sensor_geometry.auxil.p_sensor_soil = exp(-kocilai);

    # compute the fraction of sunlit leaves that can be viewed from the sensor direction (for hot spot)
    dso = sqrt( tand(CANOPY.sun_geometry.state.sza) ^ 2 +
                tand(CANOPY.sensor_geometry.state.vza) ^ 2 -
                2 * tand(CANOPY.sun_geometry.state.sza) * tand(CANOPY.sensor_geometry.state.vza) * cosd(CANOPY.sensor_geometry.state.vaa - CANOPY.sun_geometry.state.saa) );
    Σk = CANOPY.sensor_geometry.auxil.ko + CANOPY.sun_geometry.auxil.ks;
    Πk = CANOPY.sensor_geometry.auxil.ko * CANOPY.sun_geometry.auxil.ks;
    cl = CANOPY.structure.auxil.ci * CANOPY.structure.state.lai;
    α  = dso / CANOPY.structure.state.hot_spot * 2 / Σk;
    pso(x) = dso == 0 ? exp( (Σk - sqrt(Πk)) * cl * x ) : exp( Σk * cl * x + sqrt(Πk) * cl / α * (1 - exp(α * x)) );

    for i in 1:DIM_LAYER
        CANOPY.sensor_geometry.auxil.p_sun_sensor[i] = quadgk(pso, CANOPY.structure.auxil.x_bnds[i+1], CANOPY.structure.auxil.x_bnds[i]; rtol = 1e-2)[1] /
                                                       (CANOPY.structure.auxil.x_bnds[i] - CANOPY.structure.auxil.x_bnds[i+1]);
    end;

    # compute the scattering coefficients per leaf area
    for i in 1:DIM_LAYER
        j = DIM_LAYER + 1 - i;
        CANOPY.sensor_geometry.auxil.dob_leaf[:,i] .= CANOPY.sensor_geometry.auxil.dob * LEAVES[j].bio.auxil.ρ_leaf .+ CANOPY.sensor_geometry.auxil.dof * LEAVES[j].bio.auxil.τ_leaf;
        CANOPY.sensor_geometry.auxil.dof_leaf[:,i] .= CANOPY.sensor_geometry.auxil.dof * LEAVES[j].bio.auxil.ρ_leaf .+ CANOPY.sensor_geometry.auxil.dob * LEAVES[j].bio.auxil.τ_leaf;
        CANOPY.sensor_geometry.auxil.so_leaf[:,i]  .= CANOPY.sensor_geometry.auxil.sob * LEAVES[j].bio.auxil.ρ_leaf .+ CANOPY.sensor_geometry.auxil.sof * LEAVES[j].bio.auxil.τ_leaf;
    end;

    return nothing
end;

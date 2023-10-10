# This file contains functions to compute the sensor geometry of the canopy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sensor_geometry! (run per viewing zenith angle)
#
#######################################################################################################################################################################################################
"""

    sensor_geometry!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}

Update sensor geometry related auxiliary variables, given
- `config` SPAC configuration
- `can` SPAC canopy

"""
function sensor_geometry!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}
    (; Θ_AZI, Θ_INCL) = config;

    # extinction coefficients for the solar radiation
    vza = can.sensor_geometry.state.vza;
    sza = can.sun_geometry.state.sza;
    raa = can.sensor_geometry.state.vaa - can.sun_geometry.state.saa;
    for i in eachindex(Θ_INCL)
        Co = cosd(Θ_INCL[i]) * cosd(vza);
        So = sind(Θ_INCL[i]) * sind(vza);
        βo = (Co >= So ? FT(π) : acos(-Co/So));
        can.sensor_geometry.auxil.Co_incl[i] = Co;
        can.sensor_geometry.auxil.So_incl[i] = So;
        can.sensor_geometry.auxil.βo_incl[i] = βo;
        can.sensor_geometry.auxil.ko_incl[i] = 2 / FT(π) / cosd(FT(π)) * (Co * (βo - FT(π)/2) + So * sin(βo));

        # compute the scattering coefficients
        Cs = can.sun_geometry.auxil.Cs_incl[i];
        Ss = can.sun_geometry.auxil.Ss_incl[i];
        βs = can.sun_geometry.auxil.βs_incl[i];

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
        can.sensor_geometry.auxil.sb_incl[i] = (F₂ >= 0 ? F₁ : abs(F₂)) / (2 * FT(π));
        can.sensor_geometry.auxil.sf_incl[i] = (F₂ >= 0 ? F₂ : abs(F₁)) / (2 * FT(π));
    end;
    can.sensor_geometry.auxil.ko = can.structure.state.p_incl' * can.sensor_geometry.auxil.ko_incl;

    # compute the scattering weights for diffuse/direct -> sensor for backward and forward scattering
    can.sensor_geometry.auxil.dob = (can.sensor_geometry.auxil.ko + can.structure.auxil.bf) / 2;
    can.sensor_geometry.auxil.dof = (can.sensor_geometry.auxil.ko - can.structure.auxil.bf) / 2;
    can.sensor_geometry.auxil.sob = can.structure.state.p_incl' * can.sensor_geometry.auxil.sb_incl;
    can.sensor_geometry.auxil.sof = can.structure.state.p_incl' * can.sensor_geometry.auxil.sf_incl;

    # compute the fo and fo_abs matrices
    for i in eachindex(Θ_AZI)
        cos_azi_raa = cosd(Θ_AZI[i] .- (can.sensor_geometry.state.vaa - can.sun_geometry.state.saa));
        view(can.sensor_geometry.auxil.fo,:,i) .= can.sensor_geometry.auxil.Co_incl .+ can.sensor_geometry.auxil.So_incl .* cos_azi_raa;
    end;
    can.sensor_geometry.auxil.fo ./= cosd(can.sensor_geometry.state.vza);
    can.sensor_geometry.auxil.fo_abs .= abs.(can.sensor_geometry.auxil.fo);
    for i in eachindex(Θ_INCL)
        view(can.sensor_geometry.auxil.fo_cos²_incl,i,:) .= view(can.sensor_geometry.auxil.fo,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;
    can.sensor_geometry.auxil.fo_fs .= can.sun_geometry.auxil.fs .* can.sensor_geometry.auxil.fo;
    can.sensor_geometry.auxil.fo_fs_abs .= abs.(can.sensor_geometry.auxil.fo_fs);

    return nothing
end;

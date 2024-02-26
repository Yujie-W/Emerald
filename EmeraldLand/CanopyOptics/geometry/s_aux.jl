#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-25: add t_aux! method for CanopyStructureTrait-dependent variables
#
#######################################################################################################################################################################################################
s_aux!(config::SPACConfiguration{FT}, trait::CanopyStructureTrait{FT}, t_aux::CanopyStructureTDAuxil{FT}, state::SunGeometryState{FT}, s_aux::SunGeometrySDAuxil{FT}) where {FT} = (
    # if sza > 89, do nothing
    if state.sza > 89
        return nothing
    end;

    (; Θ_AZI, Θ_INCL) = config;

    # extinction coefficients for the solar radiation
    for i in eachindex(Θ_INCL)
        Cs = cosd(Θ_INCL[i]) * cosd(state.sza);
        Ss = sind(Θ_INCL[i]) * sind(state.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        s_aux.Cs_incl[i] = Cs;
        s_aux.Ss_incl[i] = Ss;
        s_aux.βs_incl[i] = βs;
        s_aux.ks_incl[i] = 2 / FT(π) / cosd(state.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;
    s_aux.ks = t_aux.p_incl' * s_aux.ks_incl;

    # compute the scattering weights for diffuse/direct -> diffuse for backward and forward scattering
    s_aux.sdb = (s_aux.ks + t_aux.bf) / 2;
    s_aux.sdf = (s_aux.ks - t_aux.bf) / 2;

    # compute the sunlit leaf fraction
    # s_aux.ps = exp.(s_aux.ks .* trait.ci * trait.lai .* can_str.t_aux.x_bnds);
    kscipai = s_aux.ks * trait.ci * (trait.lai + trait.sai);
    for i in eachindex(trait.δlai)
        ksciipai = s_aux.ks * trait.ci * (trait.δlai[i] + trait.δsai[i]);
        s_aux.p_sunlit[i] = 1 / ksciipai * (exp(kscipai * t_aux.x_bnds[i]) - exp(kscipai * t_aux.x_bnds[i+1]));
    end;

    # compute the fs and fs_abs matrices
    for i in eachindex(Θ_AZI)
        view(s_aux.fs,:,i) .= s_aux.Cs_incl .+ s_aux.Ss_incl .* cosd(Θ_AZI[i]);
    end;
    s_aux.fs ./= cosd(state.sza);
    s_aux.fs_abs .= abs.(s_aux.fs);
    mul!(s_aux.fs_abs_mean, s_aux.fs_abs', t_aux.p_incl);
    for i in eachindex(Θ_INCL)
        view(s_aux.fs_cos²_incl,i,:) .= view(s_aux.fs,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;

    return nothing
);

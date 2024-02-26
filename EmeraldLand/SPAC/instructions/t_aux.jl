t_aux!(config::SPACConfiguration{FT}, trait::CanopyStructureTrait{FT}, t_aux::CanopyStructureTDAuxil{FT}) where {FT} = (
    (; Θ_INCL_BNDS) = config;

    t_aux.x_bnds .= ([0; [sum(trait.δlai[1:i]) + sum(trait.δsai[1:i]) for i in 1:t_aux.n_layer]] ./ -(trait.lai + trait.sai));
    for i in eachindex(t_aux.p_incl)
        t_aux.p_incl[i] = lidf_cdf(trait.lidf, Θ_INCL_BNDS[i,2]) - lidf_cdf(trait.lidf, Θ_INCL_BNDS[i,1]);
    end;

    return nothing
);

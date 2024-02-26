#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-25: add t_aux! method for CanopyStructureTrait-dependent variables
#
#######################################################################################################################################################################################################
t_aux!(config::SPACConfiguration{FT}, trait::CanopyStructureTrait{FT}, t_aux::CanopyStructureTDAuxil{FT}) where {FT} = (
    (; Θ_INCL, Θ_INCL_BNDS) = config;

    # compute the LAI and SAI bounds
    t_aux.x_bnds .= ([0; [sum(trait.δlai[1:i]) + sum(trait.δsai[1:i]) for i in eachindex(trait.δlai)]] ./ -(trait.lai + trait.sai));

    # compute the probability of leaf inclination angles based on lidf
    for i in eachindex(t_aux.p_incl)
        t_aux.p_incl[i] = lidf_cdf(trait.lidf, Θ_INCL_BNDS[i,2]) - lidf_cdf(trait.lidf, Θ_INCL_BNDS[i,1]);
    end;

    # compute the weighed average of the leaf inclination angle distribution
    t_aux.bf = 0;
    for i in eachindex(Θ_INCL)
        t_aux.bf += t_aux.p_incl[i] * cosd(Θ_INCL[i]) ^ 2;;
    end;
    t_aux.ddb = (1 + t_aux.bf) / 2;
    t_aux.ddf = (1 - t_aux.bf) / 2;

    return nothing
);

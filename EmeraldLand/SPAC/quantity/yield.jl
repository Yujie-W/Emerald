#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-09: add function to compute the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p
#     2024-Jul-02: weigh ϕ_d, ϕ_f, ϕ_n, and ϕ_p by PPAR
#
#######################################################################################################################################################################################################
"""

    ϕDFNP(spac::BulkSPAC{FT}) where {FT}

Return the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p, given
- `spac` `BulkSPAC` type struct

"""
function ΦDFNP(spac::BulkSPAC{FT}) where {FT}
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    sum_ϕda::FT = 0;
    sum_ϕfa::FT = 0;
    sum_ϕna::FT = 0;
    sum_ϕpa::FT = 0;
    sum_par::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sum_ϕda += (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_d_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_d_shaded) * canopy.structure.trait.δlai[irt];
        sum_ϕfa += (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_f_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_f_shaded) * canopy.structure.trait.δlai[irt];
        sum_ϕna += (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_n_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_n_shaded) * canopy.structure.trait.δlai[irt];
        sum_ϕpa += (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_p_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_p_shaded) * canopy.structure.trait.δlai[irt];
        sum_par += (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded) * canopy.structure.trait.δlai[irt];
    end;

    return (sum_ϕda, sum_ϕfa, sum_ϕna, sum_ϕpa) ./ sum_par
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jun-06: add function to compute the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p per layer
#     2024-Jul-02: weigh ϕ_d, ϕ_f, ϕ_n, and ϕ_p by PPAR
#
#######################################################################################################################################################################################################
"""

    ΦDFNP_LAYER(spac::BulkSPAC{FT}) where {FT}

Return the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p per layer, given
- `spac` `BulkSPAC` type struct

"""
function ΦDFNP_LAYER(spac::BulkSPAC{FT}) where {FT}
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    ϕda::Vector{FT} = zeros(FT, n_layer);
    ϕfa::Vector{FT} = zeros(FT, n_layer);
    ϕna::Vector{FT} = zeros(FT, n_layer);
    ϕpa::Vector{FT} = zeros(FT, n_layer);
    ∑aa::Vector{FT} = zeros(FT, n_layer);
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        ϕda[irt] = (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_d_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_d_shaded) * canopy.structure.trait.δlai[irt];
        ϕfa[irt] = (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_f_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_f_shaded) * canopy.structure.trait.δlai[irt];
        ϕna[irt] = (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_n_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_n_shaded) * canopy.structure.trait.δlai[irt];
        ϕpa[irt] = (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit .* leaf.flux.auxil.ϕ_p_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded * leaf.flux.auxil.ϕ_p_shaded) * canopy.structure.trait.δlai[irt];
        ∑aa[irt] = (canopy.sun_geometry.s_aux.p_sunlit[irt] * mean(leaf.flux.auxil.ppar_sunlit) +
                   (1 - canopy.sun_geometry.s_aux.p_sunlit[irt]) * leaf.flux.auxil.ppar_shaded) * canopy.structure.trait.δlai[irt];
    end;

    return ϕda ./ ∑aa, ϕfa ./ ∑aa, ϕna ./ ∑aa, ϕpa ./ ∑aa
end;

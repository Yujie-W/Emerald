#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-09: add function to compute the weighted average of ϕ_f and ϕ_p
#     2024-Jul-02: weigh ϕ_f and ϕ_p by PPAR
#
#######################################################################################################################################################################################################
"""

    ΦF_ΦP(spac::BulkSPAC{FT}) where {FT}

Return the weighted average of ϕ_f and ϕ_p, given
- `spac` `BulkSPAC` type struct

"""
function ΦF_ΦP(spac::BulkSPAC{FT}) where {FT}
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    ppar_phi_f = spac.cache.cache_incl_azi_2_1;
    ppar_phi_p = spac.cache.cache_incl_azi_2_2;

    sum_ϕfa::FT = 0;
    sum_ϕpa::FT = 0;
    sum_par::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        @. ppar_phi_f = leaf.flux.auxil.ppar * leaf.photosystem.auxil.ϕ_f;
        @. ppar_phi_p = leaf.flux.auxil.ppar * leaf.photosystem.auxil.ϕ_p;
        sum_ϕfa += ppar_phi_f' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        sum_ϕpa += ppar_phi_p' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        sum_par += leaf.flux.auxil.ppar' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
    end;

    return (sum_ϕfa, sum_ϕpa) ./ sum_par
end;

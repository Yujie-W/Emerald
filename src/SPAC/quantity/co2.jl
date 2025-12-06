"""

    LEAF_PCI(spac::BulkSPAC{FT}) where {FT}

Return the weighted average of internal leaf CO₂ partial pressure, given
- `spac` `BulkSPAC` type struct

"""
function LEAF_PCI(spac::BulkSPAC{FT}) where {FT}
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    ppar_pci = spac.cache.cache_incl_azi_2_1;

    sum_pci::FT = 0;
    sum_par::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        @. ppar_pci = leaf.flux.auxil.ppar * leaf.flux.auxil.p_CO₂_i;
        sum_pci += ppar_pci' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        sum_par += leaf.flux.auxil.ppar' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
    end;

    return sum_pci / sum_par
end;

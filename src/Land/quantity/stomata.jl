"""

    BETA(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    BETA(spac::BulkSPAC{FT}) where {FT}

Return the average beta factor for
- `config` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function BETA end;

BETA(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = BETA(spac);

BETA(spac::BulkSPAC{FT}) where {FT} = (
    leaves = spac.plant.leaves;

    # compute the mean beta
    βs = 0;
    for leaf in leaves
        βs += read_β(leaf);
    end;

    return βs / length(leaves)
);


"""

    LEAF_PCI(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    LEAF_PCI(spac::BulkSPAC{FT}) where {FT}

Return the weighted average of internal leaf CO₂ partial pressure, given
- `config` `SPACConfig` type struct
- `spac` `BulkSPAC` type struct

"""
function LEAF_PCI end;

LEAF_PCI(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = LEAF_PCI(spac);

LEAF_PCI(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    ppar_pci = spac.cache.cache_incl_azi_2_1;

    sum_pci::FT = 0;
    sum_par::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        mask = leaf.flux.auxil.a_n .> 0;
        @. ppar_pci = leaf.flux.auxil.ppar * (mask * leaf.flux.auxil.p_CO₂_i);
        sum_pci += ppar_pci' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        sum_par += (leaf.flux.auxil.ppar .* mask)' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
    end;

    return sum_pci / sum_par
);

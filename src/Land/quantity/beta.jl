"""

    BETA(::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    BETA(spac::BulkSPAC{FT}) where {FT}

Return the average beta factor for
- `config` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function BETA end;

BETA(::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = BETA(spac);

BETA(spac::BulkSPAC{FT}) where {FT} = (
    leaves = spac.plant.leaves;

    # compute the mean beta
    βs = 0;
    for leaf in leaves
        βs += read_β(leaf);
    end;

    return βs / length(leaves)
);

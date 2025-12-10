"""

    CNPP(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    CNPP(spac::BulkSPAC{FT}) where {FT}

Return the canopy net primary productivity per ground area, given
- `config` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function CNPP end;

CNPP(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = CNPP(spac);

CNPP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    cnpp::FT = 0;
    for irt in eachindex(leaves)
        ilf = n_layer + 1 - irt;
        cnpp += leaves[ilf].flux.auxil.a_n_mean * canopy.structure.trait.δlai[irt];
    end;

    return cnpp
);

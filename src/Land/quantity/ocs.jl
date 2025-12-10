"""

    OCS(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    OCS(spac::BulkSPAC{FT}) where {FT}

Return the OCS flux per ground area, given
- `config` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function OCS end;

OCS(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = OCS(spac);

OCS(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    ocs::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        ocs += leaves[ilf].flux.auxil.f_ocs_mean * canopy.structure.trait.δlai[irt];
    end;

    return ocs
);

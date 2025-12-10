"""

    ET(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    ET(spac::BulkSPAC{FT}) where {FT}

Return the total evapotranspiration rate per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function ET end;

ET(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = ET(spac);

ET(spac::BulkSPAC{FT}) where {FT} = ET_SOIL(spac) + ET_VEGE(spac);


"""

    ET_SOIL(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    ET_SOIL(spac::BulkSPAC{FT}) where {FT}

Return the evaporation rate per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function ET_SOIL end;

ET_SOIL(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = ET_SOIL(spac);

ET_SOIL(spac::BulkSPAC{FT}) where {FT} = spac.soil_bulk.auxil.dndt[1,3];


"""

    ET_VEGE(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    ET_VEGE(spac::BulkSPAC{FT}) where {FT}

Return the transpiration rate per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function ET_VEGE end;

ET_VEGE(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = ET_VEGE(spac);

ET_VEGE(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute transpiration rate
    tran::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        tran += flow_out(leaves[ilf]) * canopy.structure.trait.δlai[irt];
    end;

    return tran
);

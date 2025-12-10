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


"""

    GPP(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    GPP(spac::BulkSPAC{FT}) where {FT}

Return the gross primary productivity per ground area, given
- `config` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function GPP end;

GPP(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = GPP(spac);

GPP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    gpp::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        gpp += leaves[ilf].flux.auxil.a_g_mean * canopy.structure.trait.δlai[irt];
    end;

    return gpp
);


"""

    GPP_LAYER(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    GPP_LAYER(spac::BulkSPAC{FT}) where {FT}

Return the gross primary productivity per layer, given
- `spac` `BulkSPAC` SPAC

"""
function GPP_LAYER end;

GPP_LAYER(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = GPP_LAYER(spac);

GPP_LAYER(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # get the GPP per layer
    gpps::Vector{FT} = zeros(FT, n_layer);
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        gpps[irt] = leaves[ilf].flux.auxil.a_g_mean * canopy.structure.trait.δlai[irt];
    end;

    return gpps
);


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

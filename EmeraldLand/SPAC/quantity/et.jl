#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-Jun-03: add function ET_SOIL
#
#######################################################################################################################################################################################################
"""

    ET_SOIL(spac::BulkSPAC{FT}) where {FT}

Return the evaporation rate per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ET_SOIL end;

ET_SOIL(spac::BulkSPAC{FT}) where {FT} = (
    return spac.soil_bulk.auxil.dndt[1,3]
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-08: add function to add up transpiration rate for SPAC
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    ET_VEGE(spac::BulkSPAC{FT}) where {FT}

Return the transpiration rate per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ET_VEGE end;

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

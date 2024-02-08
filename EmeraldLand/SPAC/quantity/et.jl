#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-08: add function to add up transpiration rate for SPAC
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    T_VEG(spac::BulkSPAC{FT}) where {FT}

Return the transpiration rate per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function T_VEG end;

T_VEG(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute transpiration rate
    tran::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        tran += flow_out(leaves[ilf]) * canopy.structure.state.δlai[irt];
    end;

    return tran
);

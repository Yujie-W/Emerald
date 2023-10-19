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

    # compute transpiration rate
    tran::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        tran += flow_out(leaves[i]) * canopy.structure.state.δlai[j];
    end;

    return tran
);

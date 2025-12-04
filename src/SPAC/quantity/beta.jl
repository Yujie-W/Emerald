#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-18: add function to read beta from spac
#
#######################################################################################################################################################################################################
"""

    BETA(spac::BulkSPAC{FT}) where {FT}

Return the average beta factor for
- `spac` `BulkSPAC` SPAC

"""
function BETA end;

BETA(spac::BulkSPAC{FT}) where {FT} = (
    leaves = spac.plant.leaves;

    # compute the mean beta
    βs = 0;
    for leaf in leaves
        βs += read_β(leaf);
    end;

    return βs / length(leaves)
);

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-13: add function to add up total ETR
#
#######################################################################################################################################################################################################
"""

    ΣETR(spac::BulkSPAC{FT}) where {FT}

Return the total ETR per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ΣETR end;

ΣETR(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    Σetr::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        Σetr += (canopy.sun_geometry.auxil.p_sunlit[irt] * mean(leaves[ilf].flux.auxil.etr_sunlit) +
                (1 - canopy.sun_geometry.auxil.p_sunlit[irt]) * leaves[ilf].flux.auxil.etr_shaded) * canopy.structure.state.δlai[irt];
    end;

    return Σetr
);

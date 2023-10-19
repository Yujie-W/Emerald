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

    # compute GPP
    Σetr::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        Σetr += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.etr_sunlit) +
                (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.etr_shaded) * canopy.structure.state.δlai[j];
    end;

    return Σetr
);

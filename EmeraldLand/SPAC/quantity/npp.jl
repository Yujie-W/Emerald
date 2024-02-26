#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy net primary productivity
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    CNPP(spac::BulkSPAC{FT}) where {FT}

Return the canopy net primary productivity per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function CNPP end;

CNPP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    cnpp::FT = 0;
    for irt in eachindex(leaves)
        ilf = n_layer + 1 - irt;
        cnpp += (canopy.sun_geometry.auxil.p_sunlit[irt] * mean(leaves[ilf].flux.auxil.a_n_sunlit) +
                (1 - canopy.sun_geometry.auxil.p_sunlit[irt]) * leaves[ilf].flux.auxil.a_n_shaded) * canopy.structure.trait.δlai[irt];
    end;

    return cnpp
);

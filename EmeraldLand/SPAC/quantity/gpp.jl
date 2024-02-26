#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-07: add function to add up GPP for SPAC
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    GPP(spac::BulkSPAC{FT}) where {FT}

Return the gross primary productivity per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function GPP end;

GPP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    gpp::FT = 0;
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        gpp += (canopy.sun_geometry.auxil.p_sunlit[irt] * mean(leaves[ilf].flux.auxil.a_g_sunlit) +
               (1 - canopy.sun_geometry.auxil.p_sunlit[irt]) * leaves[ilf].flux.auxil.a_g_shaded) * canopy.structure.trait.δlai[irt];
    end;

    return gpp
);

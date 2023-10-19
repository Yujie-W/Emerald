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

    # compute GPP
    gpp::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        gpp += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit) +
               (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded) * canopy.structure.state.δlai[j];
    end;

    return gpp
);

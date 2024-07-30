# This file contains functions to compute leaf photosynthesis of the entire plant

######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for BulkSPAC
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Mar-11: only compute respiration rate if solar zenith angle >= 89
#     2023-Mar-11: do nothing if LAI == 0
#     2024-Jul-25: save average a_g and a_n (to use later)
#     2024-Jul-30: compute OCS fluxes along with photosynthesis
#
#######################################################################################################################################################################################################
"""

    plant_photosynthesis!(spac::BulkSPAC{FT}) where {FT}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `BulkSPAC` type SPAC

"""
function plant_photosynthesis! end;

plant_photosynthesis!(spac::BulkSPAC{FT}) where {FT} = plant_photosynthesis!(spac, spac.plant.leaves[1]);

plant_photosynthesis!(spac::BulkSPAC{FT}, ::CanopyLayer{FT}) where {FT} = (
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    airs = spac.airs;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    n_layer = length(leaves);

    rd_only = spac.canopy.sun_geometry.state.sza > 89;
    for ilf in eachindex(leaves)
        irt = n_layer + 1 - ilf;
        leaf = leaves[ilf];
        air = airs[lindex[ilf]];
        leaf_photosynthesis!(spac.cache, leaf, air; rd_only = rd_only);

        # update the OCS flux
        leaf.flux.auxil.f_ocs .= leaf.flux.auxil.g_OCS .* air.s_aux.ps[6] ./ air.state.p_air * FT(1e6);

        # update the average rates
        leaf.flux.auxil.a_g_mean   = leaf.flux.auxil.a_g'   * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        leaf.flux.auxil.a_n_mean   = leaf.flux.auxil.a_n'   * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        leaf.flux.auxil.f_ocs_mean = leaf.flux.auxil.f_ocs' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
    end;

    return nothing
);

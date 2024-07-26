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
#
#######################################################################################################################################################################################################
"""

    plant_photosynthesis!(spac::BulkSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `BulkSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`

"""
function plant_photosynthesis! end;

plant_photosynthesis!(spac::BulkSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT} = plant_photosynthesis!(spac, mode, spac.plant.leaves[1]);

plant_photosynthesis!(spac::BulkSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}, ::CanopyLayer{FT}) where {FT} = (
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
        leaf_photosynthesis!(spac.cache, leaf, air, mode; rd_only = rd_only);  # some allocation due to broadcasting 1D to 2D array

        # update the average photosynthesis rates (a_g and a_net)
        f_sunlit = canopy.sun_geometry.s_aux.p_sunlit[irt];
        f_shaded = 1 - f_sunlit;
        leaf.flux.auxil.a_g_mean = f_sunlit * mean(view(leaf.flux.auxil.a_g,1:length(leaf.flux.auxil.a_g))) + f_shaded * leaf.flux.auxil.a_g[end];
        leaf.flux.auxil.a_n_mean = f_sunlit * mean(view(leaf.flux.auxil.a_n,1:length(leaf.flux.auxil.a_n))) + f_shaded * leaf.flux.auxil.a_n[end];
    end;

    return nothing
);

plant_photosynthesis!(spac::BulkSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}, ::Leaf{FT}) where {FT} = (
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
        leaf_photosynthesis!(spac.cache, leaf, air, mode; rd_only = rd_only);

        # update the average photosynthesis rates (a_g and a_net)
        f_sunlit = canopy.sun_geometry.s_aux.p_sunlit[irt];
        f_shaded = 1 - f_sunlit;
        leaf.flux.auxil.a_g_mean = f_sunlit * mean(leaf.flux.auxil.a_g_sunlit) + f_shaded * leaf.flux.auxil.a_g_shaded;
        leaf.flux.auxil.a_n_mean = f_sunlit * mean(leaf.flux.auxil.a_n_sunlit) + f_shaded * leaf.flux.auxil.a_n_shaded;
    end;

    return nothing
);

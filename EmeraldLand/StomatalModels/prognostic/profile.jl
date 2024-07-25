# This file contains functions to update the stomatal conductance for H₂O and CO₂

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2023-Mar-13: add function to update stomatal conductance profile based on gs and gb
#     2023-Oct-25: run nighttime stomatal conductance model if PPAR <= 1
#     2024-Jul-24: add method for matrix calculation to speed up
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance_profile!(spac::BulkSPAC{FT}) where {FT}

Compute marginal stomatal conductance change for H₂O, given
- `spac` `BulkSPAC` type struct

"""
function stomatal_conductance_profile! end;

stomatal_conductance_profile!(spac::BulkSPAC{FT}) where {FT} = (
    can_str = spac.canopy.structure;

    # if lai = 0 or roots are not connected, do nothing
    if can_str.trait.lai <= 0 || !spac.plant._root_connection
        return nothing
    end;

    airs = spac.airs;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    n_layer = length(leaves);

    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        # stomatal_conductance_profile!(leaves[ilf], airs[lindex[ilf]], can_str.auxil.ϵ_lw_layer[irt]);
        stomatal_conductance_profile!(spac.cache, leaves[ilf], airs[lindex[ilf]], can_str.auxil.ϵ_lw_layer[irt]);
    end;

    return nothing
);

stomatal_conductance_profile!(leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = (
    if leaf.flux.auxil.ppar_shaded > 1
        leaf.flux.auxil.∂g∂t_shaded = ∂g∂t(leaf, air);
        for i in eachindex(leaf.flux.auxil.∂g∂t_sunlit)
            leaf.flux.auxil.∂g∂t_sunlit[i] = ∂g∂t(leaf, air, i);
        end;
    else
        dgndt = ∂gₙ∂t(leaf, air, eff_ϵ);
        leaf.flux.auxil.∂g∂t_shaded = dgndt;
        leaf.flux.auxil.∂g∂t_sunlit .= dgndt;
    end;

    return nothing
);

stomatal_conductance_profile!(cache::SPACCache{FT}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}, eff_ϵ::FT) where {FT} = (
    if leaf.flux.auxil.ppar_shaded > 1
        leaf.flux.auxil.∂g∂t_shaded = ∂g∂t(leaf, air);
        ∂g∂t!(cache, leaf, air);
    else
        dgndt = ∂gₙ∂t(leaf, air, eff_ϵ);
        leaf.flux.auxil.∂g∂t_shaded = dgndt;
        leaf.flux.auxil.∂g∂t_sunlit .= dgndt;
    end;

    return nothing
);

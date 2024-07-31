# This file contains the functions to compute the ∂Θ∂E for WangSM optimality model

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-24: add function to update the ∂Θ∂E for sunlit leaves (matrix version)
#     2024-Jul-25: compute ∂Θ∂E for the entire leaf layer (no sunlit and shaded fractions as the sunlit and shaded part could be on the same leaf)
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E!(cache::SPACCache{FT}, sm::WangSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT}

Update the ∂Θ∂E for sunlit leaves, given
- `sm` `WangSM` type WangSM
- `leaf` `CanopyLayer` type leaf

"""
function ∂Θ∂E! end;

∂Θ∂E!(cache::SPACCache{FT}, sm::WangSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    leaf.flux.auxil.∂Θ∂E .= leaf.flux.auxil.a_n ./ max(eps(FT), (leaf.xylem.auxil.e_crit - flow_out(leaf)) / leaf.xylem.trait.area);

    return nothing
);

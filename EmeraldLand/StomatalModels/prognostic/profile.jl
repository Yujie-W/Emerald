# This file contains functions to update the stomatal conductance for H₂O and CO₂

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2023-Mar-13: add function to update stomatal conductance profile based on gs and gb
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance_profile!(spac::BulkSPAC{FT}) where {FT}

Compute marginal stomatal conductance change for H₂O, given
- `spac` `BulkSPAC` type struct

"""
function stomatal_conductance_profile! end;

stomatal_conductance_profile!(spac::BulkSPAC{FT}) where {FT} = (
    # if lai = 0 or roots are not connected, do nothing
    if spac.CANOPY.structure.state.lai == 0 || !spac._root_connection
        return nothing
    end;

    (; AIRS, LEAVES, LEAVES_INDEX) = spac;

    for i in eachindex(LEAVES_INDEX)
        stomatal_conductance_profile!(LEAVES[i], AIRS[LEAVES_INDEX[i]]);
    end;

    return nothing
);

stomatal_conductance_profile!(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    leaf.flux.auxil.∂g∂t_shaded = ∂g∂t(leaf, air);
    for i in eachindex(leaf.flux.auxil.∂g∂t_sunlit)
        leaf.flux.auxil.∂g∂t_sunlit[i] = ∂g∂t(leaf, air, i);
    end;

    return nothing
);

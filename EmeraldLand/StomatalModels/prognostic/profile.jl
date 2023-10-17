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
    if spac.canopy.structure.state.lai == 0 || !spac.plant._root_connection
        return nothing
    end;

    airs = spac.airs;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;

    for i in eachindex(leaves)
        stomatal_conductance_profile!(leaves[i], airs[lindex[i]]);
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

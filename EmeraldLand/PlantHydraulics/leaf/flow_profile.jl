# This file contains functions related to leaf flow profile

flow_in(leaf::Union{CanopyLayer{FT}, Leaf{FT}}) where {FT} = flow_in(leaf.xylem);

flow_out(leaf::Union{CanopyLayer{FT}, Leaf{FT}}) where {FT} = flow_out(leaf.xylem) + leaf.capacitor.auxil.flow;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: add function leaf_flow_profiles!
#     2023-Sep-28: account for the buffer from capacitor when running under non-steady state mode
#     2024-Feb-28: add LAI <= 0 control
#     2024-Feb-28: set leaf flow rate as the total flow rate of all area
#
#######################################################################################################################################################################################################
"""

    leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Set the flow out from each leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `BulkSPAC` type struct

"""
function leaf_flow_profiles! end;

leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = leaf_flow_profiles!(config, spac, spac.plant.leaves[1]);

leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, ::CanopyLayer{FT}) where {FT} = (
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    # run the flow profile calculation for each leaf layer only if LAI > 0
    # compute the flow rate exiting the leaf based on sunlit and shaded fractions and update it to the leaf of a BulkSPAC
    #     leaves index is from lower to upper, and thus the sunlit leaves fraction is n_layer + 1 - i
    #     airs index is also from lower to upper, but there are some layers are used by trunk so that it need to be indexed through LEAVES_INDEX

    (; ALLOW_LEAF_CONDENSATION) = config;
    airs = spac.airs;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    n_layer = length(leaves);

    g_ss = spac.cache.cache_incl_azi_1_1;

    for ilf in eachindex(leaves)
        irt = n_layer + 1 - ilf;
        leaf = leaves[ilf];

        g_ss .= 1 ./ (1 ./ leaf.flux.state.g_H₂O_s .+ 1 ./ (FT(1.35) .* leaf.flux.auxil.g_CO₂_b));
        g = g_ss' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        d = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000) - airs[lindex[ilf]].s_aux.ps[3];
        ALLOW_LEAF_CONDENSATION ? nothing : d = max(d, 0);
        f = g * d / airs[lindex[ilf]].state.p_air * leaf.xylem.trait.area;

        # set_flow_out!
        set_flow_profile!(leaf.xylem, f - leaf.capacitor.auxil.flow);
    end;

    return nothing
);

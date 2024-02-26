# This file contains functions related to leaf flow profile

flow_in(leaf::Leaf{FT}) where {FT} = flow_in(leaf.xylem);

flow_out(leaf::Leaf{FT}) where {FT} = flow_out(leaf.xylem) + leaf.capacitor.auxil.flow;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: add function leaf_flow_profiles!
#     2023-Sep-28: account for the buffer from capacitor when running under non-steady state mode
#
#######################################################################################################################################################################################################
"""

    leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Set the flow out from each leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `BulkSPAC` type struct

"""
function leaf_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # compute the flow rate exiting the leaf based on sunlit and shaded fractions and update it to the leaf of a BulkSPAC
    #     leaves index is from lower to upper, and thus the sunlit leaves fraction is n_layer + 1 - i
    #     airs index is also from lower to upper, but there are some layers are used by trunk so that it need to be indexed through LEAVES_INDEX

    (; ALLOW_LEAF_CONDENSATION) = config;
    airs = spac.airs;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    n_layer = length(leaves);

    for i in eachindex(leaves)
        leaf = leaves[i];
        f_sl = canopy.sun_geometry.s_aux.p_sunlit[n_layer + 1 - i];

        g_sh = 1 / (1 /leaf.flux.state.g_H₂O_s_shaded + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
        g_sl = 0;
        for j in eachindex(leaf.flux.state.g_H₂O_s_sunlit)
            g_sl += 1 / (1 / leaf.flux.state.g_H₂O_s_sunlit[j] + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
        end;
        g_sl /= length(leaf.flux.state.g_H₂O_s_sunlit);

        g = g_sh * (1 - f_sl) + g_sl * f_sl;
        d = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.auxil.p_leaf * 1000000) - airs[lindex[i]].s_aux.ps[3];
        ALLOW_LEAF_CONDENSATION ? nothing : d = max(d, 0);
        f = g * d / airs[lindex[i]].state.p_air;

        # set_flow_out!
        set_flow_profile!(leaf.xylem, f - leaf.capacitor.auxil.flow);
    end;

    return nothing
end;

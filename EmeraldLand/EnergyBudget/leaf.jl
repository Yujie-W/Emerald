# This file contains function to calculate energy budgets of the leaf

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_energy_flows!(spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the leaf, given
- `spac` `BulkSPAC` type SPAC

"""
function leaf_energy_flows!(spac::BulkSPAC{FT}) where {FT}
    branches = spac.plant.branches;
    airs = spac.airs;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    canopy = spac.canopy;
    # the total energy change of the leaf is the sum of
    #     the energy of the flow from the stem
    #     the energy of the flow from the air
    #     the sensible heat from the air
    #     the net radiation energy from shortwave and longwave radiation

    N = length(branches)
    for i in eachindex(branches)
        stem = branches[i];
        leaf = leaves[i];
        air = airs[lindex[i]];

        # if flow in is positive, then energy flow is positive for leaf
        f_i = flow_in(leaf);
        if f_i >= 0
            leaf.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * stem.energy.auxil.t;
        else
            leaf.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * leaf.energy.auxil.t;
        end;

        # if flow out is positive, then energy flow is negative for leaf
        # if f_o is positive, then the leaf is losing water to the air
        # note here that CP_L_MOL is included in the latent_heat_vapor TD function
        # so here we only need to calculate the heat mass flow from water vapor in gas phase
        f_o = flow_out(leaf);
        leaf.energy.auxil.∂e∂t -= f_o * M_H₂O(FT) * latent_heat_vapor(leaf.energy.auxil.t);
        if f_o >= 0
            leaf.energy.auxil.∂e∂t -= f_o * CP_V_MOL(FT) * leaf.energy.auxil.t;
        else
            leaf.energy.auxil.∂e∂t -= f_o * CP_V_MOL(FT) * air.auxil.t;
        end;

        # add the sensible heat flux from the leaf to air (to total leaf area)
        g_be = FT(1.4) * FT(0.135) * sqrt(air.auxil.wind / (FT(0.72) * leaf.bio.state.width));
        leaf.energy.auxil.∂e∂t -= 2 * g_be * CP_D_MOL(FT) * (leaf.energy.auxil.t - air.auxil.t) * leaf.xylem.state.area;

        # add the net radiation energy to the leaf (to total leaf area)
        leaf.energy.auxil.∂e∂t += (canopy.sun_geometry.auxil.r_net_sw[N+1-i] + canopy.structure.auxil.r_net_lw_leaf[N+1-i]) / canopy.structure.state.δlai[N+1-i] * leaf.xylem.state.area;
    end;

    return nothing
end;

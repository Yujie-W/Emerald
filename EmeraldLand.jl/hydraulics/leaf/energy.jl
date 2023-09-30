# This file contains function to calculate energy budgets of the leaf

heat_capacitance(leaf::Leaf{FT}) where {FT} = heat_capacitance(
    # leaf heat capaciatance is the sum of the heat capaciatance of the water stored in the leaf and the heat capaciatance of the leaf itself
    # here convert lma from g cm⁻² to kg m⁻² with the factor 10
    return leaf.capacitor.state.v_storage * CP_L_MOL(FT) + leaf.xylem.state.cp * leaf.xylem.state.area * leaf.bio.state.lma * 10
);

heat_capacitance(leaf::Leaves2D{FT}) where {FT} = heat_capacitance(leaf.NS);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_energy_flows!(spac::MultiLayerSPAC{FT}) where {FT}

Calculate the energy flows of the leaf, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function leaf_energy_flows!(spac::MultiLayerSPAC{FT}) where {FT}
    (; AIR, BRANCHES, LEAVES, LEAVES_INDEX) = spac;
    # the total energy change of the leaf is the difference between
    #     the energy of the flow from the stem
    #     the energy of the flow to the air

    for i in eachindex(BRANCHES)
        stem = BRANCHES[i];
        leaf = LEAVES[i];
        air = AIR[LEAVES_INDEX[i]];

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
            leaf.energy.auxil.∂e∂t -= f_o * CP_V_MOL(FT) * air.t;
        end;
    end;

    return nothing
end;

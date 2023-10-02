# This file contains functions to calculate the energy budgets of the junctions

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the junction
#
#######################################################################################################################################################################################################
"""

    junction_energy_flows!(spac::MultiLayerSPAC{FT}) where {FT}

Calculate the energy flows of the junction, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function junction_energy_flows!(spac::MultiLayerSPAC{FT}) where {FT}
    (; JUNCTION, ROOTS, TRUNK) = spac;

    # The total energy change of the junction is difference between
    #     the energy of the flow from the roots
    #     the energy of the flow to the trunk
    # if the flow out of root is positive, then the energy flow is positive
    for root in ROOTS
        f_o = flow_out(root);
        if f_o >= 0
            JUNCTION.auxil.∂e∂t += f_o * CP_L_MOL(FT) * root.energy.auxil.t;
        else
            JUNCTION.auxil.∂e∂t += f_o * CP_L_MOL(FT) * JUNCTION.auxil.t;
        end;
    end;

    # if the flow into the trunk is positive, then the energy flow is negative
    f_i = flow_in(TRUNK);
    if f_i >= 0
        JUNCTION.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * JUNCTION.auxil.t;
    else
        JUNCTION.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * TRUNK.energy.auxil.t;
    end;

    return nothing
end;

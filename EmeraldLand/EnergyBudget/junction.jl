# This file contains functions to calculate the energy budgets of the junctions

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the junction
#
#######################################################################################################################################################################################################
"""

    junction_energy_flows!(spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the junction, given
- `spac` `BulkSPAC` type SPAC

"""
function junction_energy_flows!(spac::BulkSPAC{FT}) where {FT}
    junction = spac.plant.junction;
    roots = spac.plant.roots;
    trunk = spac.plant.trunk;

    # The total energy change of the junction is difference between
    #     the energy of the flow from the roots
    #     the energy of the flow to the trunk
    # if the flow out of root is positive, then the energy flow is positive
    for root in roots
        f_o = flow_out(root);
        if f_o >= 0
            junction.auxil.∂e∂t += f_o * CP_L_MOL(FT) * root.energy.s_aux.t;
        else
            junction.auxil.∂e∂t += f_o * CP_L_MOL(FT) * junction.s_aux.t;
        end;
    end;

    # if the flow into the trunk is positive, then the energy flow is negative
    f_i = flow_in(trunk);
    if f_i >= 0
        junction.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * junction.s_aux.t;
    else
        junction.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * trunk.energy.s_aux.t;
    end;

    return nothing
end;

# This file contains function to calculate energy budgets of the stem

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the stem (trunk and branches)
#     2023-Oct-18: add net radiation to branch energy budget
#
#######################################################################################################################################################################################################
"""

    stem_energy_flows!(spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the trunk and branches, given
- `spac` `BulkSPAC` type SPAC

"""
function stem_energy_flows!(spac::BulkSPAC{FT}) where {FT}
    branches = spac.plant.branches;
    canopy = spac.canopy;
    junction = spac.plant.junction;
    leaves = spac.plant.leaves;
    sbulk = spac.soil_bulk;
    trunk = spac.plant.trunk;
    n_layer = length(leaves);

    # for the trunk, the total energy is the differentce of
    #     energy from the junction
    #     energy to the branches
    # if flow in is positive, then energy flow is positive
    f_i = flow_in(trunk);
    if f_i >= 0
        trunk.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * junction.auxil.t;
    else
        trunk.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * trunk.energy.auxil.t;
    end;

    # for the branches, the total energy is the difference of
    #     energy from the trunk
    #     energy to the leaves
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        stem = branches[ilf];
        leaf = leaves[ilf];

        # if flow in is positive, then energy flow is positive for branches but negative for trunk
        f_i = flow_in(stem);
        if f_i >= 0
            stem.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * trunk.energy.auxil.t;
            trunk.energy.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * trunk.energy.auxil.t;
        else
            stem.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * stem.energy.auxil.t;
            trunk.energy.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * stem.energy.auxil.t;
        end;

        # if flow out is positive, then energy flow is negative for branchesß
        f_o = flow_out(stem);
        if f_o >= 0
            stem.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * stem.energy.auxil.t;
        else
            stem.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * leaf.energy.auxil.t;
        end;

        # add the net radiation energy to the leaf (to total leaf area)
        stem.energy.auxil.∂e∂t += (canopy.sun_geometry.auxil.r_net_sw_stem[irt] + canopy.structure.auxil.r_net_lw_stem[irt]) * sbulk.state.area;
    end;

    return nothing
end;

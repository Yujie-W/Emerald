# This file contains function to calculate energy budgets of the stem

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the stem (trunk and branches)
#
#######################################################################################################################################################################################################
"""

    stem_energy_flows!(spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the trunk and branches, given
- `spac` `BulkSPAC` type SPAC

"""
function stem_energy_flows!(spac::BulkSPAC{FT}) where {FT}
    (; BRANCHES, JUNCTION, LEAVES, TRUNK) = spac;

    # for the trunk, the total energy is the differentce of
    #     energy from the junction
    #     energy to the branches
    # if flow in is positive, then energy flow is positive
    f_i = flow_in(TRUNK);
    if f_i >= 0
        TRUNK.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * JUNCTION.auxil.t;
    else
        TRUNK.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * TRUNK.energy.auxil.t;
    end;

    # for the branches, the total energy is the difference of
    #     energy from the trunk
    #     energy to the leaves
    for i in eachindex(BRANCHES)
        stem = BRANCHES[i];
        leaf = LEAVES[i];

        # if flow in is positive, then energy flow is positive for branches but negative for trunk
        f_i = flow_in(stem);
        if f_i >= 0
            stem.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * TRUNK.energy.auxil.t;
            TRUNK.energy.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * TRUNK.energy.auxil.t;
        else
            stem.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * stem.energy.auxil.t;
            TRUNK.energy.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * stem.energy.auxil.t;
        end;

        # if flow out is positive, then energy flow is negative for branchesß
        f_o = flow_out(stem);
        if f_o >= 0
            stem.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * stem.energy.auxil.t;
        else
            stem.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * leaf.energy.auxil.t;
        end;
    end;

    return nothing
end;

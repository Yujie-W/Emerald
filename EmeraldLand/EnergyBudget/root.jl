# This file contains function to calculate energy budgets of the root

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add function to calculate the energy flow of the root
#     2023-Oct-06: add function to calculate the energy flow between root and soil
#
#######################################################################################################################################################################################################
"""

    root_energy_flows!(spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the root, given
- `spac` `BulkSPAC` type SPAC

"""
function root_energy_flows!(spac::BulkSPAC{FT}) where {FT}
    junction = spac.plant.junction;
    rindx = spac.plant.roots_index;
    roots = spac.plant.roots;
    sbulk = spac.soil_bulk;
    soils = spac.soils;

    # compute the energy flow per layer
    # the energy flow is computed as the difference between
    #     energy of the water flow from soil
    #     energy of the water flow to the root-trunk junction
    for i in eachindex(roots)
        root = roots[i];
        soil = soils[rindx[i]];

        # if the flow into the root is positive, then the energy flow is positive
        f_i = flow_in(root);
        if f_i >= 0
            root.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * soil.auxil.t;
            soil.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * soil.auxil.t / sbulk.state.area;
        else
            root.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * root.energy.auxil.t;
            soil.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * root.energy.auxil.t / sbulk.state.area;
        end;

        # if the flow into the junction is positive, then the energy flow is negative
        f_o = flow_out(root);
        if f_o >= 0
            root.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * root.energy.auxil.t;
        else
            root.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * junction.auxil.t;
        end;
    end;

    return nothing
end;

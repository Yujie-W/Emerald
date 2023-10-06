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

    root_energy_flows!(spac::MultiLayerSPAC{FT}) where {FT}

Calculate the energy flows of the root, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function root_energy_flows!(spac::MultiLayerSPAC{FT}) where {FT}
    (; JUNCTION, ROOTS, ROOTS_INDEX, SOIL_BULK, SOILS) = spac;

    # compute the energy flow per layer
    # the energy flow is computed as the difference between
    #     energy of the water flow from soil
    #     energy of the water flow to the root-trunk junction
    for i in eachindex(ROOTS)
        root = ROOTS[i];
        soil = SOILS[ROOTS_INDEX[i]];

        # if the flow into the root is positive, then the energy flow is positive
        f_i = flow_in(root);
        if f_i >= 0
            root.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * soil.auxil.t;
            soil.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * soil.auxil.t / SOIL_BULK.state.area;
        else
            root.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * root.energy.auxil.t;
            soil.auxil.∂e∂t -= f_i * CP_L_MOL(FT) * root.energy.auxil.t / SOIL_BULK.state.area;
        end;

        # if the flow into the junction is positive, then the energy flow is negative
        f_o = flow_out(root);
        if f_o >= 0
            root.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * root.energy.auxil.t;
        else
            root.energy.auxil.∂e∂t -= f_o * CP_L_MOL(FT) * JUNCTION.auxil.t;
        end;
    end;

    return nothing
end;

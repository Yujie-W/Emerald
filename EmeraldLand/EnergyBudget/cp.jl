# This file contains functions to calculate the energy budgets of the xylem (root and stem)

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-29: add function heat_capacitance
#     2023-Oct-06: add method for SoilLayer
#     2023-Oct-07: account for the runoff water in the heat capacitance of the top soil
#     2023-Oct-09: add method for AirLayer
#
#######################################################################################################################################################################################################
"""

    heat_capacitance(organ)

Return the heat capacitance of the organ (xylem, root, stem, and leaf)

"""
function heat_capacitance end;

heat_capacitance(soil::SoilLayer{FT}; runoff::FT = FT(0)) where {FT} = (
    cp_gas = (soil.state.ns[3] * CP_V_MOL(FT) + (soil.state.ns[1] + soil.state.ns[2] + soil.state.ns[4] + soil.state.ns[5]) * CP_D_MOL(FT)) / soil.auxil.δz;

    # runoff in mol m⁻² s⁻¹, convert it to kg and then

    return soil.state.ρ * soil.state.cp + soil.state.θ * ρ_H₂O(FT) * CP_L(FT) + cp_gas + runoff * CP_L_MOL(FT) / soil.auxil.δz
);

heat_capacitance(xylem::XylemHydraulics{FT}) where {FT} = (
    # The heat capaciatance of the xylem is the sum of the heat capaciatance of the water stored in the xylem and the heat capaciatance of the xylem itself
    return sum(xylem.state.v_storage) * CP_L_MOL(FT) + xylem.state.cp * xylem.state.area * xylem.state.l
);

heat_capacitance(root::Root{FT}) where {FT} = heat_capacitance(root.xylem);

heat_capacitance(junc::JunctionCapacitor{FT}) where {FT} = (
    # the heat capaciatance of the junction is that of all water stored in the junction
    return junc.state.v_storage * CP_L_MOL(FT)
);

heat_capacitance(root::Stem{FT}) where {FT} = heat_capacitance(root.xylem);

heat_capacitance(leaf::Leaf{FT}) where {FT} = heat_capacitance(
    # leaf heat capaciatance is the sum of the heat capaciatance of the water stored in the leaf and the heat capaciatance of the leaf itself
    # here convert lma from g cm⁻² to kg m⁻² with the factor 10
    return leaf.capacitor.state.v_storage * CP_L_MOL(FT) + leaf.xylem.state.cp * leaf.xylem.state.area * leaf.bio.state.lma * 10
);

heat_capacitance(air::AirLayer{FT}) where {FT} = (air.state.p_air - air.auxil.ps[3]) * air.auxil.δz / GAS_R(FT) * CP_D_MOL(FT) + air.state.ns[3] * CP_V_MOL(FT);

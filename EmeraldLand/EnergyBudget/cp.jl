# This file contains functions to calculate the energy budgets of the xylem (root and stem)

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-29: add function heat_capacitance
#
#######################################################################################################################################################################################################
"""

    heat_capacitance(organ)

Return the heat capacitance of the organ (xylem, root, stem, and leaf)

"""
function heat_capacitance end;

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

heat_capacitance(leaf::Leaves2D{FT}) where {FT} = heat_capacitance(leaf.NS);

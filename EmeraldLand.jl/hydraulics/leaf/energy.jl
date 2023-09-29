# This file contains function to calculate energy budgets of the leaf

heat_capacitance(leaf::Leaf{FT}) where {FT} = heat_capacitance(
    # leaf heat capaciatance is the sum of the heat capaciatance of the water stored in the leaf and the heat capaciatance of the leaf itself
    # here convert lma from g cm⁻² to kg m⁻² with the factor 10
    return leaf.capacitor.state.v_storage * CP_L_MOL(FT) + leaf.xylem.state.cp * leaf.xylem.state.area * leaf.bio.state.lma * 10
);

heat_capacitance(leaf::Leaves2D{FT}) where {FT} = heat_capacitance(leaf.NS.);

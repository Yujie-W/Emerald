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

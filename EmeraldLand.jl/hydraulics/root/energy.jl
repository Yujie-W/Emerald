# This file contains function to calculate energy budgets of the root

heat_capacitance(root::Root{FT}) where {FT} = heat_capacitance(root.xylem);

heat_capacitance(root::Root2{FT}) where {FT} = heat_capacitance(root.NS.xylem);

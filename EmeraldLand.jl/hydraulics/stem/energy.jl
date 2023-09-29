# This file contains function to calculate energy budgets of the stem

heat_capacitance(root::Stem{FT}) where {FT} = heat_capacitance(root.xylem);

heat_capacitance(root::Stem2{FT}) where {FT} = heat_capacitance(root.NS.xylem);

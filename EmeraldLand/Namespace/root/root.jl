# This file contains the root struct as a combination of the xylem and energy structs
# Note here that the energy struct needs to be initialized with the heat capacity of root and water, water content, and temperature...

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add Root struct with energy and xylem fields
#     2023-Sep-23: add constructor for Root struct and initialize the energy state of the root
#     2023-Sep-23: add field rhizosphere
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Root{FT}
    "Root energy struct"
    energy::XylemEnergy{FT} = XylemEnergy{FT}()
    "Rhizosphere struct"
    rhizosphere::Rhizosphere{FT} = Rhizosphere{FT}()
    "Root xylem struct"
    xylem::XylemHydraulics{FT}
end;


"""

    Root(config::SPACConfiguration{FT}) where {FT}

Return the root struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Root(config::SPACConfiguration{FT}) where {FT} = (
    root = Root{FT}(xylem = XylemHydraulics(config));
    initialize_energy_storage!(root);

    return root
);

initialize_energy_storage!(root::Root{FT}) where {FT} = (
    root.xylem.state.v_storage .= (root.xylem.state.v_max * root.xylem.state.area * root.xylem.state.l) / length(root.xylem.state.v_storage);
    root.energy.auxil.cp = sum(root.xylem.state.v_storage) * CP_L_MOL(FT) + (root.xylem.state.cp * root.xylem.state.area * root.xylem.state. l);
    root.energy.state.Σe = root.energy.auxil.cp * root.energy.auxil.t;

    return nothing
);

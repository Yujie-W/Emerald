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
    energy::XylemEnergy{FT}
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
    r_energy = XylemEnergy{FT}();
    r_xylem = XylemHydraulics(config);

    # now update the energy state of the root before returning the root struct
    r_energy.auxil.cp = sum(r_xylem.state.v_storage) * CP_L_MOL(FT) + (r_xylem.state.cp * r_xylem.state.area * r_xylem.state. l);
    r_energy.state.energy = r_energy.auxil.cp * r_energy.auxil.t;

    return Root{FT}(
                energy = r_energy,
                xylem  = r_xylem
    )
);

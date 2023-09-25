# This file contains stem state and auxil structs

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-25: add new Stem struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save stem energy and hydraulics variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Stem{FT}
    "Stem energy struct"
    energy::XylemEnergy{FT}
    "Stem hydraulics struct"
    xylem::XylemHydraulics{FT}
end


"""

    Stem(config::SPACConfiguration{FT}) where {FT}

Return the stem struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Stem(config::SPACConfiguration{FT}) where {FT} = (
    s_energy = XylemEnergy{FT}();
    s_xylem = XylemHydraulics(config);

    # now update the energy state of the stem before returning the stem struct
    s_energy.auxil.cp = sum(s_xylem.state.v_storage) * CP_L_MOL(FT) + (s_xylem.state.cp * s_xylem.state.area * s_xylem.state. l);
    s_energy.state.energy = s_energy.auxil.cp * s_energy.auxil.t;

    return Stem{FT}(
                energy = s_energy,
                xylem  = s_xylem
    )
);

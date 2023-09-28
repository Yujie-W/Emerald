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
    energy::XylemEnergy{FT} = XylemEnergy{FT}()
    "Stem hydraulics struct"
    xylem::XylemHydraulics{FT}
end


"""

    Stem(config::SPACConfiguration{FT}) where {FT}

Return the stem struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Stem(config::SPACConfiguration{FT}) where {FT} = (
    stem = Stem{FT}(xylem = XylemHydraulics(config));
    initialize_energy_storage!(stem);

    return stem
);

initialize_energy_storage!(stem::Stem{FT}) where {FT} = (
    stem.xylem.state.v_storage .= (stem.xylem.state.v_max * stem.xylem.state.area * stem.xylem.state.l) / length(stem.xylem.state.v_storage);
    stem.energy.auxil.cp = sum(stem.xylem.state.v_storage) * CP_L_MOL(FT) + (stem.xylem.state.cp * stem.xylem.state.area * stem.xylem.state. l);
    stem.energy.state.energy = stem.energy.auxil.cp * stem.energy.auxil.t;

    return nothing
);

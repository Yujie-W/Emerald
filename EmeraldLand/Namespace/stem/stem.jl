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
Stem(config::SPACConfiguration{FT}) where {FT} = Stem{FT}(xylem = XylemHydraulics(config));

# This file contains the rhizosphere information

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add RhizosphereState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save rhizosphere state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RhizosphereState{FT}
    "Maximum rhizosphere conductance per root area `[mol s⁻¹ m⁻² MPa⁻¹]`"
    k_max::FT = 1e20
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add RhizosphereAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save rhizosphere auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RhizosphereAuxil{FT}
    "Rhizosphere pressure `[mol s⁻¹ MPa⁻¹]`"
    p_rhizo::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add Rhizosphere
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save rhizosphere parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Rhizosphere{FT}
    "Rhizosphere state"
    state::RhizosphereState{FT} = RhizosphereState{FT}()
    "Rhizosphere auxilary variables"
    auxil::RhizosphereAuxil{FT} = RhizosphereAuxil{FT}()
end;

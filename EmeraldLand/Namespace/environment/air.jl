# This file contains the state and auxiliary variable for the air layer

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct AirLayerState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayerState{FT}
    "Moles of CH₄, CO₂, H₂O, N₂, and O₂ `[mol]`"
    ns::Vector{FT} = FT[0, 0, 0, 0, 0]
    "Atmospheric pressure `[Pa]`"
    p_air::FT = P_ATM(FT)
    "The lower and upper boundary of the air layer `[m]`"
    zs::Vector{FT} = FT[0, 1]
    "Total energy within the air layer `[J m⁻²]`"
    Σe::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct AirLayerAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayerAuxil{FT}
    "CO₂ concentration `[ppm]`"
    f_CO₂::FT = 400
    "Partial pressures of CH₄, CO₂, H₂O, N₂, and O₂ `[Pa]`"
    ps::Vector{FT} = FT[0, 0, 1500, 0, 0]
    "Temperature `[K]`"
    t::FT = T₂₅(FT)
    "Wind speed `[m s⁻¹]`"
    wind::FT = 1
    "The mean location of the air layer"
    z::FT = 0.5
    "The thickness of the air layer"
    δz::FT = 1
    "Marginal increase in total energy `[J m⁻² s⁻¹]`"
    ∂e∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: refactor the structure to contain only state and auxiliary struct fields
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    "State variables"
    state::AirLayerState{FT} = AirLayerState{FT}()
    "Auxiliary variables"
    auxil::AirLayerAuxil{FT} = AirLayerAuxil{FT}()
end;

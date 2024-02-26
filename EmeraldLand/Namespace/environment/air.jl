# This file contains the state and auxiliary variable for the air layer

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct AirLayerTrait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer trait variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayerTrait{FT}
    "The lower and upper boundary of the air layer `[m]`"
    zs::Vector{FT} = FT[0, 1]
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct AirLayerTDAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer trait derived auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayerTDAuxil{FT}
    "The thickness of the air layer"
    δz::FT = 1
    "The mean location of the air layer"
    z::FT = 0.5
end;

t_aux!(trait::AirLayerTrait{FT}, t_aux::AirLayerTDAuxil{FT}) where {FT} = (
    t_aux.δz = trait.zs[2] - trait.zs[1];
    t_aux.z = (trait.zs[1] + trait.zs[2]) / 2;

    return nothing
);


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
    "Total energy within the air layer `[J m⁻²]`"
    Σe::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct AirLayerSDAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer state derived auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayerSDAuxil{FT}
    "CO₂ concentration `[ppm]`"
    f_CO₂::FT = 400
    "Partial pressures of CH₄, CO₂, H₂O, N₂, and O₂ `[Pa]`"
    ps::Vector{FT} = FT[0, 40, 1500, 0, 0]
    "Temperature `[K]`"
    t::FT = T₂₅(FT)
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
    "Wind speed `[m s⁻¹]`"
    wind::FT = 1
    "Marginal increase in total energy `[J m⁻² s⁻¹]`"
    ∂e∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: refactor the structure to contain only state and auxiliary struct fields
#     2024-Feb-26: add field trait, t_aux, and s_aux
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    "Trait variables"
    trait::AirLayerTrait{FT} = AirLayerTrait{FT}()
    "State variables"
    state::AirLayerState{FT} = AirLayerState{FT}()
    "Trait-derived auxiliary variables"
    t_aux::AirLayerTDAuxil{FT} = AirLayerTDAuxil{FT}()
    "State-derived auxiliary variables"
    s_aux::AirLayerSDAuxil{FT} = AirLayerSDAuxil{FT}()
    "Auxiliary variables"
    auxil::AirLayerAuxil{FT} = AirLayerAuxil{FT}()
end;

t_aux!(air::AirLayer{FT}) where {FT} = t_aux!(air.trait, air.t_aux);

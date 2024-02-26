# This file contains the soil layer state and auxiliary variables

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct SoilLayerTrait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer trait variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayerTrait{FT}
    "Specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    cp::FT = 760
    "Soil moisture retention curve"
    vc::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")
    "Soil thermal conductivity `[W m⁻¹ K⁻¹]`"
    λ_soil::FT = 3
    "Dry soil density `[kg m⁻³]`"
    ρ::FT = 2650

    # Geometry information
    "Depth boundaries `[m]`"
    zs::Vector{FT} = FT[0,-1]
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct SoilLayerTDAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer trait dependent auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayerTDAuxil{FT}
    "Mean depth `[m]`"
    z::FT = 0.5
    "Layer thickness `[m]`"
    δz::FT = 1
end;

t_aux!(trait::SoilLayerTrait{FT}, t_aux::SoilLayerTDAuxil{FT}) where {FT} = (
    t_aux.z = abs(trait.zs[1] + trait.zs[2]) / 2;
    t_aux.δz = trait.zs[1] - trait.zs[2];

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilLayerState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayerState{FT}
    "Moles of soil trace gasses (CH₄, CO₂, H₂O, N₂, O₂) `[mol]`"
    ns::Vector{FT} = zeros(FT,5)
    "Soil water content"
    θ::FT = 0.3
    "Total stored energy per volume `[J m⁻³]`"
    Σe::FT = 0
end;


Base.@kwdef mutable struct SoilLayerSDAuxil{FT}
    "Combined specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    cp::FT = 0
    "Temperature `[K]`"
    t::FT = T₂₅(FT)
    "Soil hydraulic conductance per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    k::FT = 0
    "Relative soil diffusive coefficient per area based on air fraction (distance accounted for already)"
    kd::FT = 0
    "Relative soil diffusive coefficient for water vapor (distance accounted for already)"
    kv::FT = 0
    "Combined soil thermal conductance `[W m⁻² K⁻¹]`"
    λ_soil_water::FT = 0
    "Matric potential `[MPa]`"
    ψ::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilLayerAuxil
#     2023-Oct-05: add field kv
#     2023-Oct-05: add field n_con
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayerAuxil{FT}
    "Moles of condensated water vapor `[mol]`"
    n_con::FT = 0
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase in trace gas moles `[mol s⁻¹]`"
    ∂n∂t::Vector{FT} = zeros(FT,5)
    "Marginal increase in soil water content `[s⁻¹]`"
    ∂θ∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilLayer
#     2024-Feb-26: add field trait, t_aux, and s_aux
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer state and auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayer{FT}
    "Trait variables"
    trait::SoilLayerTrait{FT} = SoilLayerTrait{FT}()
    "State variables"
    state::SoilLayerState{FT} = SoilLayerState{FT}()
    "Trait dependent auxiliary variables"
    t_aux::SoilLayerTDAuxil{FT} = SoilLayerTDAuxil{FT}()
    "State dependent auxiliary variables"
    s_aux::SoilLayerSDAuxil{FT} = SoilLayerSDAuxil{FT}()
    "Auxiliary variables"
    auxil::SoilLayerAuxil{FT} = SoilLayerAuxil{FT}()
end;

t_aux!(layer::SoilLayer{FT}) where {FT} = t_aux!(layer.trait, layer.t_aux);

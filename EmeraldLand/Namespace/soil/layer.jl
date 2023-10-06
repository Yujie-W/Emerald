# This file contains the soil layer state and auxiliary variables

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
    "Specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    cp::FT = 760
    "Moles of soil trace gasses (CH₄, CO₂, H₂O, N₂, O₂) `[mol]`"
    ns::Vector{FT} = zeros(FT,5)
    "Soil moisture retention curve"
    vc::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")
    "Soil thermal conductivity `[W m⁻¹ K⁻¹]`"
    λ_soil::FT = 3
    "Dry soil density `[kg m⁻³]`"
    ρ::FT = 2650
    "Soil water content"
    θ::FT = vc.Θ_SAT
    "Total stored energy per volume `[J m⁻³]`"
    Σe::FT = 0

    # Geometry information
    "Depth boundaries `[m]`"
    zs::Vector{FT} = FT[0,-1]
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilLayerAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayerAuxil{FT}
    "Combined specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    cp::FT = 0
    "Soil hydraulic conductance per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    k::FT = 0
    "Relative soil diffusive coefficient per area based on air fraction (distance accounted for already)"
    kd::FT = 0
    "Temperature `[K]`"
    t::FT = T₂₅(FT)
    "Mean depth `[m]`"
    z::FT = 0.5
    "Layer thickness `[m]`"
    δz::FT = 1
    "Combined soil thermal conductance `[W m⁻² K⁻¹]`"
    λ_soil_water::FT = 0
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase in trace gas moles `[mol s⁻¹]`"
    ∂n∂t::Vector{FT} = zeros(FT,5)
    "Marginal increase in soil water content `[s⁻¹]`"
    ∂θ∂t::FT = 0
    "Matric potential `[MPa]`"
    ψ::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilLayer
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil layer state and auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayer{FT}
    "State variables"
    state::SoilLayerState{FT} = SoilLayerState{FT}()
    "Auxiliary variables"
    auxil::SoilLayerAuxil{FT} = SoilLayerAuxil{FT}()
end;

initialize_energy_storage!(layer::SoilLayer{FT}) where {FT} = (
    layer.state.Σe = layer.state.ρ * layer.state.cp * layer.auxil.t + layer.state.θ * CP_L(FT) * ρ_H₂O(FT) * layer.auxil.t;

    return nothing;
);

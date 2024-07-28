# This file contains the state and auxil variables related to stomtal conductance (carbon and water fluxes)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-25: add CanopyLayerFluxState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf flux state variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyLayerFluxState{FT}
    "Stomatal conductance to water vapor for sunlit and shaded (end element) leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s::Vector{FT}
end;

CanopyLayerFluxState(config::SPACConfiguration{FT}) where {FT} = CanopyLayerFluxState{FT}(0.01 .* ones(FT, config.DIM_PPAR_BINS+1));


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-25: add CanopyLayerFluxAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf flux auxiliary variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyLayerFluxAuxil{FT}
    # stomtal conductance
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Total leaf diffusive conductance to CO₂ for sunlit and shaded (end element) leaves `[mol m⁻² s⁻¹]`"
    g_CO₂::Vector{FT}
    "Marginal increase of conductance per time for sunlit and shaded (end element) leaves `[mol m⁻² s⁻²]`"
    ∂g∂t::Vector{FT}
    "Marginal increase in A per increase in transpiration rate for sunlit and shaded (end element) leaves `[μmol m⁻² s⁻¹]`"
    ∂A∂E::Vector{FT}
    "Marginal increase in Θ per increase in transpiration rate"
    ∂Θ∂E::FT = 0

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure for sunlit and shaded (end element) leaves `[Pa]`"
    p_CO₂_i::Vector{FT}
    "Leaf surface CO₂ partial pressure for sunlit and shaded (end element) leaves `[Pa]`"
    p_CO₂_s::Vector{FT}

    # Photosynthesis
    "Gross photosynthetic rate for sunlit and shaded (end element) leaves `[μmol m⁻² s⁻¹]`"
    a_g::Vector{FT}
    "Average gross photosynthetic rate [μmol m⁻² s⁻¹]"
    a_g_mean::FT = 0
    "Average net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_n_mean::FT = 0
    "Net photosynthetic rate for sunlit and shaded (end element) leaves `[μmol m⁻² s⁻¹]`"
    a_n::Vector{FT}

    # Integrators
    "Integrator for transpiration out"
    ∫∂w∂t_out::FT = 0

    # ppar from canopy radiation
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit and shaded (end element) leaves `[μmol m⁻² s⁻¹]`"
    ppar::Vector{FT}

    # used for nocturnal stomatal conductance
    "Memory PPAR `[μmol m⁻² s⁻¹]`"
    ppar_mem::FT = 1000

    # used for empirical model
    "Beta of the empirical models (NaN for optimality models)"
    β::FT = NaN
end;

CanopyLayerFluxAuxil(config::SPACConfiguration{FT}) where {FT} = CanopyLayerFluxAuxil{FT}(
    g_CO₂   = zeros(FT, config.DIM_PPAR_BINS+1),
    ∂g∂t    = zeros(FT, config.DIM_PPAR_BINS+1),
    ∂A∂E    = zeros(FT, config.DIM_PPAR_BINS+1),
    p_CO₂_i = zeros(FT, config.DIM_PPAR_BINS+1),
    p_CO₂_s = zeros(FT, config.DIM_PPAR_BINS+1),
    a_g     = zeros(FT, config.DIM_PPAR_BINS+1),
    a_n     = zeros(FT, config.DIM_PPAR_BINS+1),
    ppar    = zeros(FT, config.DIM_PPAR_BINS+1),
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add sturct LeafFlux
#     2024-Feb-26: add fiel trait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf flux variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyLayerFlux{FT}
    "Trait variables"
    trait::LeafFluxTrait{FT} = LeafFluxTrait{FT}()
    "Leaf flux state variables"
    state::CanopyLayerFluxState{FT}
    "Leaf flux auxiliary variables"
    auxil::CanopyLayerFluxAuxil{FT}
end;

CanopyLayerFlux(config::SPACConfiguration{FT}) where {FT} = CanopyLayerFlux{FT}(state = CanopyLayerFluxState(config), auxil = CanopyLayerFluxAuxil(config));

# This file contains the state and auxil variables related to stomtal conductance (carbon and water fluxes)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add sturct LeafFluxTrait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf flux trait variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafFluxTrait{FT}
    "Stomtal model"
    stomatal_model::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    g_limits::Vector{FT} = FT[1e-4, 0.3]
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add sturct LeafFluxState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf flux state variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafFluxState{FT}
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::FT = 0.01
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add sturct LeafFluxAuxil
#     2023-Oct-24: add fields ϕ_f1_* and ϕ_f2_*
#     2024-Jul-24: add field ∂A∂E and ∂Θ∂E to compute the dgdt using matrix calculation (much faster)
#     2024-Aug-05: remove the shaded and sunlit distinctions as this is now meant for leaf level only
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf flux auxiliary variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafFluxAuxil{FT}
    # stomtal conductance
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::FT = 0
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::FT = 0
    "Marginal increase in A per increase in transpiration rate"
    ∂A∂E::FT = 0
    "Marginal increase in Θ per increase in transpiration rate"
    ∂Θ∂E::FT = 0

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i::FT = 0
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_CO₂_s::FT = 0

    # Integrators
    "Integrator for transpiration out"
    ∫∂w∂t_out::FT = 0

    # ppar from canopy radiation
    "Absorbed photosynthetically active radiation `[μmol m⁻² s⁻¹]`"
    apar::FT = 0
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::FT = 0

    # used for nocturnal stomatal conductance
    "Memory PPAR `[μmol m⁻² s⁻¹]`"
    ppar_mem::FT = 1000

    # used for empirical model
    "Beta of the empirical models (NaN for optimality models)"
    β::FT = NaN
end;


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
Base.@kwdef mutable struct LeafFlux{FT}
    "Trait variables"
    trait::LeafFluxTrait{FT} = LeafFluxTrait{FT}()
    "Leaf flux state variables"
    state::LeafFluxState{FT} = LeafFluxState{FT}()
    "Leaf flux auxiliary variables"
    auxil::LeafFluxAuxil{FT} = LeafFluxAuxil{FT}()
end;

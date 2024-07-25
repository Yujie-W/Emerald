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
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT = 0.01
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT}
end;

LeafFluxState(config::SPACConfiguration{FT}) where {FT} = LeafFluxState{FT}(g_H₂O_s_sunlit = 0.01 .* ones(FT, config.DIM_INCL, config.DIM_AZI));


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add sturct LeafFluxAuxil
#     2023-Oct-24: add fields ϕ_f1_* and ϕ_f2_*
#     2024-Jul-24: add field ∂A∂E_sunlit and ∂Θ∂E to compute the dgdt using matrix calculation (much faster)
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
    "Total leaf diffusive conductance to CO₂ for shaded leaves `[mol m⁻² s⁻¹]`"
    g_CO₂_shaded::FT = 0
    "Total leaf diffusive conductance to CO₂ for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_CO₂_sunlit::Matrix{FT}
    "Marginal increase of conductance per time for shaded leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_shaded::FT = 0
    "Marginal increase of conductance per timefor sunlit leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_sunlit::Matrix{FT}
    "Marginal increase in A per increase in transpiration rate for shaded leaves"
    ∂A∂E_shaded::FT = 0
    "Marginal increase in A per increase in transpiration rate for sunlit leaves"
    ∂A∂E_sunlit::Matrix{FT}
    "Marginal increase in Θ per increase in transpiration rate"
    ∂Θ∂E::FT = 0

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure for shaded leaves `[Pa]`"
    p_CO₂_i_shaded::FT = 0
    "Leaf internal CO₂ partial pressure for sunlit leaves `[Pa]`"
    p_CO₂_i_sunlit::Matrix{FT}
    "Leaf surface CO₂ partial pressure for shaded leaves `[Pa]`"
    p_CO₂_s_shaded::FT = 0
    "Leaf surface CO₂ partial pressure for sunlit leaves `[Pa]`"
    p_CO₂_s_sunlit::Matrix{FT}

    # Photosynthesis
    "Average gross photosynthetic rate [μmol m⁻² s⁻¹]"
    a_g::FT = 0
    "Gross photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_g_shaded::FT = 0
    "Gross photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_g_sunlit::Matrix{FT}
    "Average net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_n::FT = 0
    "Net photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_n_shaded::FT = 0
    "Net photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_n_sunlit::Matrix{FT}
    "Actual electron transport for shaded leaves `[μmol m⁻² s⁻¹]`"
    etr_shaded::FT = 0
    "Actual electron transport for sunlit leaves `[μmol m⁻² s⁻¹]`"
    etr_sunlit::Matrix{FT}
    "Heat dissipation quantum yield for shaded leaves `[-]`"
    ϕ_d_shaded::FT = 0
    "Heat dissipation quantum yield for sunlit leaves `[-]`"
    ϕ_d_sunlit::Matrix{FT}
    "Fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f_shaded::FT = 0
    "Fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f_sunlit::Matrix{FT}
    "Non-photochemical quenching quantum yield for shaded leaves `[-]`"
    ϕ_n_shaded::FT = 0
    "Non-photochemical quenching quantum yield for sunlit leaves `[-]`"
    ϕ_n_sunlit::Matrix{FT}
    "Photochemical quantum yield for shaded leaves `[-]`"
    ϕ_p_shaded::FT = 0
    "Photochemical quantum yield for sunlit leaves `[-]`"
    ϕ_p_sunlit::Matrix{FT}

    # Fluorescence yeilds of two photosystems
    "PSI fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f1_shaded::FT = 0
    "PSI fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f1_sunlit::Matrix{FT}
    "PSII fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f2_shaded::FT = 0
    "PSII fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f2_sunlit::Matrix{FT}

    # Integrators
    "Integrator for transpiration out"
    ∫∂w∂t_out::FT = 0

    # ppar from canopy radiation
    "Absorbed photosynthetically active radiation for shaded leaves `[μmol m⁻² s⁻¹]`"
    apar_shaded::FT = 0
    "Absorbed photosynthetically active radiation for sunlit leaves `[μmol m⁻² s⁻¹]`"
    apar_sunlit::Matrix{FT}
    "Absorbed photosynthetically active radiation used for photosynthesis for shaded leaves `[μmol m⁻² s⁻¹]`"
    ppar_shaded::FT = 0
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit leaves `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Matrix{FT}

    # used for nocturnal stomatal conductance
    "Memory PPAR `[μmol m⁻² s⁻¹]`"
    ppar_mem::FT = 1000

    # used for empirical model
    "Beta of the empirical models (NaN for optimality models)"
    β::FT = NaN
end;

LeafFluxAuxil(config::SPACConfiguration{FT}) where {FT} = LeafFluxAuxil{FT}(
            g_CO₂_sunlit   = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ∂g∂t_sunlit    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ∂A∂E_sunlit    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            p_CO₂_i_sunlit = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            p_CO₂_s_sunlit = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            a_g_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            a_n_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            etr_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ϕ_d_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ϕ_f_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ϕ_n_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ϕ_p_sunlit     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ϕ_f1_sunlit    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ϕ_f2_sunlit    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            apar_sunlit    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            ppar_sunlit    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
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
Base.@kwdef mutable struct LeafFlux{FT}
    "Trait variables"
    trait::LeafFluxTrait{FT} = LeafFluxTrait{FT}()
    "Leaf flux state variables"
    state::LeafFluxState{FT}
    "Leaf flux auxiliary variables"
    auxil::LeafFluxAuxil{FT}
end;

LeafFlux(config::SPACConfiguration{FT}) where {FT} = LeafFlux{FT}(state = LeafFluxState(config), auxil = LeafFluxAuxil(config));

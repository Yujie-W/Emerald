#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-27: add new structure for leaves with 2D Matrix of parameters for sunlit partitioning and point value for shaded partitioning
#     2022-Jun-27: make BIO HyperspectralLeafBiophysics only
#     2022-Jun-27: add sunlit and shaded ppar to struct (remove the ppar in canopy radiation)
#     2022-Jun-28: add a_gross, a_net, and ϕ_f for sunlit and shaded leaves
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#     2022-Jul-12: add fields: ∂g∂t_shaded and ∂g∂t_sunlit
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Jul-19: add dimension control to struct
#     2022-Jul-28: move field _t to PSM
#     2022-Nov-18: use Union type for SM
#     2023-Mar-02: set minimum G to 1e-4 instead of 1e-2
#     2023-Apr-13: move field APAR_CAR to SPACConfiguration
#     2023-Jun-13: add field: etr_shaded, etr_sunlit
#     2023-Jun-16: remove fields DIM_*
#     2023-Sep-07: add water flow integrators
#     2023-Sep-09: add fields ϕ_x_shaded and ϕ_x_sunlit
#     2023-Sep-11: set minimum G to 0 instead of 1e-4
#     2023-Sep-18: use HyperLeafBio instead of HyperspectralLeafBiophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters for a single canopy layer. This structure is meant for canopy level research and canopy radiative transfer scheme with sunlit and shaded partitioning as well as leaf
    angular distribution.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaves2D{FT}
    # Constants
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-3, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "New leaf struct, will replace Leaves2D in the next major refactor"
    NS::Leaf{FT}
    # "[`HyperLeafBio`](@ref) type leaf biophysical parameters"
    # BIO::HyperLeafBio{FT}
    # "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    # HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Absorbed photosynthetically active radiation used for photosynthesis for shaded leaves `[μmol m⁻² s⁻¹]`"
    ppar_shaded::FT = 200
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit leaves `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Matrix{FT} =

    # Prognostic variables (used for ∂y∂t)
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT = 0.01
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT}
    "Marginal increase of conductance per time for shaded leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_shaded::FT = 0
    "Marginal increase of conductance per timefor sunlit leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_sunlit::Matrix{FT}

    # Diagnostic variables
    "Gross photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_gross_shaded::FT = 0
    "Gross photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_gross_sunlit::Matrix{FT}
    "Net photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_net_shaded::FT = 0
    "Net photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_net_sunlit::Matrix{FT}
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
    "Integrator for transpiration out"
    ∫∂w∂t_out = 0

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::FT = 0
    "Total leaf diffusive conductance to CO₂ for shaded leaves `[mol m⁻² s⁻¹]`"
    _g_CO₂_shaded::FT = 0
    "Total leaf diffusive conductance to CO₂ for sunlit leaves `[mol m⁻² s⁻¹]`"
    _g_CO₂_sunlit::Matrix{FT}
    "Leaf internal CO₂ partial pressure for shaded leaves `[Pa]`"
    _p_CO₂_i_shaded::FT = 0
    "Leaf internal CO₂ partial pressure for sunlit leaves `[Pa]`"
    _p_CO₂_i_sunlit::Matrix{FT}
    "Leaf surface CO₂ partial pressure for shaded leaves `[Pa]`"
    _p_CO₂_s_shaded::FT = 0
    "Leaf surface CO₂ partial pressure for sunlit leaves `[Pa]`"
    _p_CO₂_s_sunlit::Matrix{FT}
end

Leaves2D(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_AZI, DIM_INCL) = config;

    return Leaves2D{FT}(
                NS              = Leaf(config),
                ppar_sunlit     = 1000 .* ones(FT, DIM_INCL, DIM_AZI),
                g_H₂O_s_sunlit  = FT(0.01) .* ones(FT, DIM_INCL, DIM_AZI),
                ∂g∂t_sunlit     = zeros(FT, DIM_INCL, DIM_AZI),
                a_gross_sunlit  = zeros(FT, DIM_INCL, DIM_AZI),
                a_net_sunlit    = zeros(FT, DIM_INCL, DIM_AZI),
                etr_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_d_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_f_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_n_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_p_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                _g_CO₂_sunlit   = zeros(FT, DIM_INCL, DIM_AZI),
                _p_CO₂_i_sunlit = zeros(FT, DIM_INCL, DIM_AZI),
                _p_CO₂_s_sunlit = zeros(FT, DIM_INCL, DIM_AZI),
    )
);

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
    # Embedded structures
    "New leaf struct, will replace Leaves2D in the next major refactor"
    NS::Leaf{FT}
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

#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Root structure
#     2022-Jul-15: add fields e, ∂e∂t
#     2023-Sep-07: add water flow integrators
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Root2{FT<:AbstractFloat}
    # Embedded structures
    "New root struct, will replace Root2 in the next major refactor"
    NS::Root{FT}
    # "[`RootHydraulics`](@ref) type root hydraulic system"
    # HS::RootHydraulics{FT} = RootHydraulics{FT}()

    # Prognostic variables (not used for ∂y∂t)
    # "Current temperature `[K]`"
    # t::FT = T₂₅(FT)

    # Prognostic variables (used for ∂y∂t)
    # "Total stored energy in water `[J]`" # TODO: add wood storage as well
    # e::FT = sum(HS.v_storage) * CP_L_MOL(FT) * t
    # "Marginal increase in energy `[W]`"
    # ∂e∂t::FT = 0
    "Integrator for soil water extration `[mol m⁻² s⁻¹]`"
    ∫∂w∂t_in::FT = 0
    "Integrator for root water extration `[mol m⁻² s⁻¹]`"
    ∫∂w∂t_out::FT = 0

    # Cache variable
    "Whether root is connected to soil"
    _isconnected::Bool = true
end

Root2(config::SPACConfiguration{FT}) where {FT} = Root2{FT}(
            NS = Root(config)
);

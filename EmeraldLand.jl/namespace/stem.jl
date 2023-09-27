#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Stem structure
#     2022-Jul-15: add fields e, ∂e∂t
#     2022-Jul-19: use kwdef for the constructor
#     2023-Sep-07: add water flow integrators
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save stem parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Stem2{FT<:AbstractFloat}
    # Embedded structures
    "New stem struct, will replace Stem2 in the next major refactor"
    NS::Stem{FT}
    # "[`StemHydraulics`](@ref) type stem hydraulic system"
    # HS::StemHydraulics{FT} = StemHydraulics{FT}()

    # Prognostic variables (not used for ∂y∂t)
    # "Current temperature"
    # t::FT = T₂₅(FT)

    # Prognostic variables (used for ∂y∂t)
    # "Total stored energy in water `[J]`" # TODO: add wood storage as well
    # e::FT = T₂₅(FT) * sum(HS.v_storage) * CP_L_MOL(FT)
    # "Marginal increase in energy `[W]`"
    # ∂e∂t::FT = 0
    "Integrator for downstream sap water extration `[mol m⁻² s⁻¹]`"
    ∫∂w∂t_in::FT = 0
    "Integrator for stem water extration `[mol m⁻² s⁻¹]`"
    ∫∂w∂t_out::FT = 0
end

Stem2(config::SPACConfiguration{FT}) where {FT} = Stem2{FT}(
            NS = Stem(config)
);

#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Root structure
#     2022-Jul-15: add fields e, ∂e∂t
#     2023-Apr-13: add DIM_XYLEM to struct type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Root{FT,DIMS}
    # Embedded structures
    "[`RootHydraulics`](@ref) type root hydraulic system"
    HS::RootHydraulics{FT,DIMS} = RootHydraulics{FT,DIMS}()

    # Prognostic variables (not used for ∂y∂t)
    "Current temperature `[K]`"
    t::FT = T₂₅()

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy in water `[J]`" # TODO: add wood storage as well
    e::FT = sum(HS.v_storage) * CP_L_MOL(FT) * t
    "Marginal increase in energy `[W]`"
    ∂e∂t::FT = 0

    # Cache variable
    "Whether root is connected to soil"
    _isconnected::Bool = true
end

Root(config::SPACConfiguration{FT,DIMS}) where {FT,DIMS} = Root{FT,DIMS}();

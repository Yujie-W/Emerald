#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Root structure
#     2022-Jul-15: add fields e, ∂e∂t
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Root{FT}
    # Embedded structures
    "[`RootHydraulics`](@ref) type root hydraulic system"
    HS::RootHydraulics{FT,DIM_XYLEM} = RootHydraulics{FT,DIM_XYLEM}()

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

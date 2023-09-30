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

    # Cache variable
    "Whether root is connected to soil"
    _isconnected::Bool = true
end

Root2(config::SPACConfiguration{FT}) where {FT} = Root2{FT}(
            NS = Root(config)
);

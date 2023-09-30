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
end

Stem2(config::SPACConfiguration{FT}) where {FT} = Stem2{FT}(
            NS = Stem(config)
);

# This file contains the structs for canopy structural parameters

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructureState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural state variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructureState{FT}
    # canopy structure
    "Hot spot parameter"
    hot_spot::FT = 0.05
    "Leaf inclination angle distribution function algorithm"
    lidf::Union{BetaLIDF{FT}, VerhoefLIDF{FT}} = VerhoefLIDF{FT}()
    "Inclination angle distribution"
    p_incl::Vector{FT}

    # Leaf area index
    "Leaf area index"
    lai::FT
    "Leaf area index distribution"
    δlai::Vector{FT}

    # Clumping index
    "Clumping structure a"
    Ω_A::FT = 1
    "Clumping structure b"
    Ω_B::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructureAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural auxiliary variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructureAuxil{FT}
    "Weighted sum of cos²(inclination)"
    bf::FT = 0
    "Clumping index"
    ci::FT = 1
    "Canopy level boundary locations"
    x_bnds::Vector{FT}
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructure{FT}
    "State variables"
    state::CanopyStructureState{FT}
    "Auxiliary variables"
    auxil::CanopyStructureAuxil{FT}
end;

CanopyStructure(config::SPACConfiguration{FT}) where {FT} = (
    lai = 3;
    δlai = 3 .* ones(FT, config.DIM_LAYER) ./ config.DIM_LAYER;
    x_bnds = ([0; [sum(δlai[1:i]) for i in 1:config.DIM_LAYER]] ./ -lai);
    p_incl = ones(FT, config.DIM_INCL) ./ config.DIM_INCL;

    return CanopyStructure{FT}(
                state = CanopyStructureState{FT}(p_incl = p_incl, lai = lai, δlai = δlai),
                auxil = CanopyStructureAuxil{FT}(x_bnds = x_bnds),
    )
);

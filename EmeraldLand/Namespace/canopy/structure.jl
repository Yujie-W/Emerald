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
    "Clumping index"
    ci::FT = 1
    "Canopy level boundary locations"
    x_bnds::Vector{FT}

    # canopy scattering coefficients
    "Weighted sum of cos²(inclination)"
    bf::FT = 0
    "Backward diffuse->diffuse scatter weight"
    ddb::FT = 0
    "Forward diffuse->diffuse scatter weight"
    ddf::FT = 0

    # Longwave radiation (sun and sensor independent)
    "Effective emissivity for different layers"
    ϵ_lw_layer::Vector{FT}
    "Reflectance for longwave radiation at each canopy layer"
    ρ_lw_layer::Vector{FT}
    "Tranmittance for longwave radiation at each canopy layer"
    τ_lw_layer::Vector{FT}
    "Effective reflectance for longwave radiation"
    ρ_lw::Vector{FT}
    "Effective tranmittance for longwave radiation"
    τ_lw::Vector{FT}
end;

CanopyStructureAuxil(config::SPACConfiguration{FT}) where {FT} = CanopyStructureAuxil{FT}(
            x_bnds     = zeros(FT, config.DIM_LAYER + 1),
            ϵ_lw_layer = zeros(FT, config.DIM_LAYER),
            ρ_lw_layer = zeros(FT, config.DIM_LAYER),
            τ_lw_layer = zeros(FT, config.DIM_LAYER),
            ρ_lw       = zeros(FT, config.DIM_LAYER + 1),
            τ_lw       = zeros(FT, config.DIM_LAYER),
);


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
    p_incl = ones(FT, config.DIM_INCL) ./ config.DIM_INCL;

    cs_auxil = CanopyStructureAuxil(config);
    cs_auxil.x_bnds .= ([0; [sum(δlai[1:i]) for i in 1:config.DIM_LAYER]] ./ -lai);

    return CanopyStructure{FT}(
                state = CanopyStructureState{FT}(p_incl = p_incl, lai = lai, δlai = δlai),
                auxil = cs_auxil,
    )
);

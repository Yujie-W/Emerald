# This file contains the structs for canopy structural parameters

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructureState
#     2023-Oct-18: add fields sai and δsai
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

    # Leaf area index
    "Leaf area index"
    lai::FT
    "Leaf area index distribution"
    δlai::Vector{FT}

    # Stem area index
    "Stem area index"
    sai::FT
    "Stem area index distribution"
    δsai::Vector{FT}

    # Clumping index
    "Clumping structure a"
    Ω_A::FT = 1
    "Clumping structure b"
    Ω_B::FT = 0
end;

sync_state!(state_from::CanopyStructureState{FT}, state_to::CanopyStructureState{FT}) where {FT} = (
    state_to.hot_spot = state_from.hot_spot;
    sync_state!(state_from.lidf, state_to.lidf);
    state_to.lai = state_from.lai;
    state_to.δlai .= state_from.δlai;
    state_to.sai = state_from.sai;
    state_to.δsai .= state_from.δsai;
    state_to.Ω_A = state_from.Ω_A;
    state_to.Ω_B = state_from.Ω_B;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructureAuxil
#     2023-Oct-18: add fields ddb_stem, ddf_stem, lw_layer_leaf, lw_layer_stem, r_net_lw_leaf, r_net_lw_stem
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
    "Inclination angle distribution"
    p_incl::Vector{FT}
    "Canopy level boundary locations"
    x_bnds::Vector{FT}

    # canopy scattering coefficients
    "Weighted sum of cos²(inclination)"
    bf::FT = 0
    "Backward diffuse->diffuse scatter weight"
    ddb::FT = 0
    "Forward diffuse->diffuse scatter weight"
    ddf::FT = 0

    # Scattering coefficients per leaf area
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of leaf"
    ddb_leaf::Matrix{FT}
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of leaf"
    ddf_leaf::Matrix{FT}
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of stem"
    ddb_stem::Matrix{FT}
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of stem"
    ddf_stem::Matrix{FT}

    # Reflectance and tranmittance per canopy layer (no denominator correction made yet)
    "Reflectance for diffuse->diffuse at each canopy layer"
    ρ_dd_layer::Matrix{FT}
    "Tranmittance for diffuse->diffuse at each canopy layer"
    τ_dd_layer::Matrix{FT}

    # Effective reflectance and tranmittance per canopy layer (including the denominator correction)
    "Effective reflectance for diffuse->diffuse"
    ρ_dd::Matrix{FT}
    "Effective tranmittance for diffuse->diffuse"
    τ_dd::Matrix{FT}

    # Longwave radiation coefficients (sun and sensor independent)
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

    # Longwave radiation flux
    "Longwave energy flux from leaves and stem (one side) `[W m⁻²]`"
    lw_layer::Vector{FT}
    "Longwave energy flux from leaves (one side) `[W m⁻²]`"
    lw_layer_leaf::Vector{FT}
    "Longwave energy flux from stem (one side) `[W m⁻²]`"
    lw_layer_stem::Vector{FT}
    "Downwelling longwave energy flux `[W m⁻²]`"
    lwꜜ::Vector{FT}
    "Upwelling longwave energy flux `[W m⁻²]`"
    lwꜛ::Vector{FT}
    "Downwelling longwave energy flux `[W m⁻²]`"
    emitꜜ::Vector{FT}
    "Upwelling longwave energy flux `[W m⁻²]`"
    emitꜛ::Vector{FT}

    # Net longwave radiation flux
    "Net longwave energy absorption per leaf area `[W m⁻²]`"
    r_net_lw_leaf::Vector{FT}
    "Net longwave energy absorption per stem area `[W m⁻²]`"
    r_net_lw_stem::Vector{FT}
end;

CanopyStructureAuxil(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = CanopyStructureAuxil{FT}(
            p_incl        = zeros(FT, config.DIM_INCL),
            x_bnds        = zeros(FT, n_layer + 1),
            ddb_leaf      = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddf_leaf      = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddb_stem      = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddf_stem      = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ρ_dd_layer    = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            τ_dd_layer    = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ρ_dd          = zeros(FT, length(config.SPECTRA.Λ), n_layer + 1),
            τ_dd          = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ϵ_lw_layer    = zeros(FT, n_layer),
            ρ_lw_layer    = zeros(FT, n_layer),
            τ_lw_layer    = zeros(FT, n_layer),
            ρ_lw          = zeros(FT, n_layer + 1),
            τ_lw          = zeros(FT, n_layer),
            lw_layer      = zeros(FT, n_layer),
            lw_layer_leaf = zeros(FT, n_layer),
            lw_layer_stem = zeros(FT, n_layer),
            lwꜜ           = zeros(FT, n_layer + 1),
            lwꜛ           = zeros(FT, n_layer + 1),
            emitꜜ         = zeros(FT, n_layer),
            emitꜛ         = zeros(FT, n_layer + 1),
            r_net_lw_leaf = zeros(FT, n_layer),
            r_net_lw_stem = zeros(FT, n_layer),
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

CanopyStructure(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = (
    lai = 3;
    δlai = 3 .* ones(FT, n_layer) ./ n_layer;
    sai = 0.5;
    δsai = 0.5 .* ones(FT, n_layer) ./ n_layer;

    cs_auxil = CanopyStructureAuxil(config, n_layer);
    cs_auxil.x_bnds .= ([0; [sum(δlai[1:i]) + sum(δsai[1:i]) for i in 1:n_layer]] ./ -(lai + sai));
    cs_auxil.p_incl = ones(FT, config.DIM_INCL) ./ config.DIM_INCL;

    return CanopyStructure{FT}(
                state = CanopyStructureState{FT}(lai = lai, δlai = δlai, sai = sai, δsai = δsai),
                auxil = cs_auxil,
    )
);

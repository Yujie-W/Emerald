# This file contains the structs for canopy structural parameters

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructureState
#     2023-Oct-18: add fields sai and δsai
#     2024-Feb-27: move hot_spot parameter to SPACConfiguration
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: set ci to be a ClumpingIndex type struct to account for angular dependency
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural trait variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructureTrait{FT}
    # canopy structure
    "Leaf inclination angle distribution function algorithm"
    lidf::Union{BetaLIDF{FT}, VerhoefLIDF{FT}} = VerhoefLIDF{FT}()
    "Stem inclination angle distribution function algorithm"
    sidf::Union{BetaLIDF{FT}, VerhoefLIDF{FT}} = VerhoefLIDF{FT}()

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

    # Clumping index of the canopy
    "Clumping index"
    ci::ClumpingIndex{FT} = ClumpingIndex{FT}()
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-25: add struct CanopyStructureTDAuxil
#     2024-Mar-01: add field kd for diffuse radiation extinction coefficient (purely unabsorbed, not after reabsorption)
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: add field ci_diffuse
#     2024-Oct-16: change kd_leaf and kd_stem to be 90 elements long (per degree of a isotropic radiation)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural trait-dependent auxiliary variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructureTDAuxil{FT}
    # Angles
    "Inclination angle distribution of leaves"
    p_incl_leaf::Vector{FT}
    "Inclination angle distribution of stems"
    p_incl_stem::Vector{FT}
    "Canopy level boundary locations"
    x_bnds::Vector{FT}

    # Clumping index for diffuse radiation
    "Clumping index for diffuse radiation"
    ci_diffuse::FT = 1

    # diffuse radiation extinction coefficients
    "Diffuse radiation extinction coefficient (leaf; weighed with leaf angle and clumping index; per direction)"
    kd_leaf::Vector{FT} = zeros(FT, 90)
    "Diffuse radiation extinction coefficient (stem; weighed with leaf angle and clumping index; per direction)"
    kd_stem::Vector{FT} = zeros(FT, 90)

    # canopy scattering coefficients (leaf)
    "Backward diffuse->diffuse scatter weight (leaf)"
    ddb_leaf::FT = 0
    "Forward diffuse->diffuse scatter weight (leaf)"
    ddf_leaf::FT = 0

    # canopy scattering coefficients (stem)
    "Backward diffuse->diffuse scatter weight (stem)"
    ddb_stem::FT = 0
    "Forward diffuse->diffuse scatter weight (stem)"
    ddf_stem::FT = 0
end;

CanopyStructureTDAuxil(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = CanopyStructureTDAuxil{FT}(
            p_incl_leaf = ones(FT, config.DIM_INCL) ./ config.DIM_INCL,
            p_incl_stem = ones(FT, config.DIM_INCL) ./ config.DIM_INCL,
            x_bnds = zeros(FT, n_layer + 1),
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructureAuxil
#     2023-Oct-18: add fields ddb_stem, ddf_stem, lw_layer_leaf, lw_layer_stem, r_net_lw_leaf, r_net_lw_stem
#     2024-Oct-16: add field τ_dd_isotropic
#     2024-Oct-16: add fields ρ_leaf_eff and τ_leaf_eff
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural auxiliary variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructureAuxil{FT}
    # Effective leaf reflectance and transmittance for solar radiation
    "Effective leaf reflectance after accounting for the CI effect"
    ρ_leaf_eff::Matrix{FT}
    "Effective leaf transmittance after accounting for the CI  effect"
    τ_leaf_eff::Matrix{FT}

    # Scattering coefficients per leaf area
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of leaf"
    ddb_leaf::Matrix{FT}
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of leaf"
    ddf_leaf::Matrix{FT}
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of stem"
    ddb_stem::Matrix{FT}
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins of stem"
    ddf_stem::Matrix{FT}

    # Weight tranmittance
    "The transmittance chance that isotropic radiation will not reach any leaf surface"
    τ_dd_isotropic::Vector{FT}

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
            ρ_leaf_eff     = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            τ_leaf_eff     = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddb_leaf       = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddf_leaf       = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddb_stem       = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ddf_stem       = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            τ_dd_isotropic = zeros(FT, n_layer),
            ρ_dd_layer     = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            τ_dd_layer     = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ρ_dd           = zeros(FT, length(config.SPECTRA.Λ), n_layer + 1),
            τ_dd           = zeros(FT, length(config.SPECTRA.Λ), n_layer),
            ϵ_lw_layer     = zeros(FT, n_layer),
            ρ_lw_layer     = zeros(FT, n_layer),
            τ_lw_layer     = zeros(FT, n_layer),
            ρ_lw           = zeros(FT, n_layer + 1),
            τ_lw           = zeros(FT, n_layer),
            lw_layer       = zeros(FT, n_layer),
            lw_layer_leaf  = zeros(FT, n_layer),
            lw_layer_stem  = zeros(FT, n_layer),
            lwꜜ            = zeros(FT, n_layer + 1),
            lwꜛ            = zeros(FT, n_layer + 1),
            emitꜜ          = zeros(FT, n_layer),
            emitꜛ          = zeros(FT, n_layer + 1),
            r_net_lw_leaf  = zeros(FT, n_layer),
            r_net_lw_stem  = zeros(FT, n_layer),
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-09: add struct CanopyStructure
#     2024-Feb-25: add field trait, t_aux (remove state)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores canopy structural variables.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyStructure{FT}
    "Trait variables that need to be presribed from GriddingMachine"
    trait::CanopyStructureTrait{FT}
    "Trait-dependent variables"
    t_aux::CanopyStructureTDAuxil{FT}
    "Auxiliary variables"
    auxil::CanopyStructureAuxil{FT}
end;

CanopyStructure(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = (
    lai = 3;
    δlai = 3 .* ones(FT, n_layer) ./ n_layer;
    sai = 0.5;
    δsai = 0.5 .* ones(FT, n_layer) ./ n_layer;

    trait = CanopyStructureTrait{FT}(lai = lai, δlai = δlai, sai = sai, δsai = δsai);
    t_aux = CanopyStructureTDAuxil(config, n_layer);
    auxil = CanopyStructureAuxil(config, n_layer);
    t_aux.x_bnds .= ([0; [sum(δlai[1:i]) + sum(δsai[1:i]) for i in 1:n_layer]] ./ -(lai + sai));
    t_aux.p_incl_leaf = ones(FT, config.DIM_INCL) ./ config.DIM_INCL;

    return CanopyStructure{FT}(trait = trait, t_aux = t_aux, auxil = auxil)
);

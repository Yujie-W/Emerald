"""
Struct that contains the trait and state variables.
"""
Base.@kwdef mutable struct LeafTraitsStates{FT<:AbstractFloat}
    # leaf pigments
    "Anthocyanin content `[μg cm⁻²]`"
    ant::FT = 0
    "Senescent material (brown pigments) fraction `[-]`"
    brown::FT = 0
    "Chlorophyll a and b content `[μg cm⁻²]`"
    cab::FT = 40
    "Carotenoid content `[μg cm⁻²]`"
    car::FT = 40 / 7
    "Carbon-based constituents in lma `[g cm⁻²]`"
    cbc::FT = 0
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
    "Leaf mesophyll structural parameter that describes the number of thin layers with a leaf"
    meso_n::FT = 1.4
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT = 0

    # leaf pigments related states
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::FT = 0
    "Fraction of carotenoid absorption into SIF `[-]`"
    ϕ_car::FT = 1
    "Fraction of carotenoid absorption into PPAR `[-]`"
    ϕ_car_ppar::FT = 1

    # leaf width
    "leaf width `[m]`"
    width::FT = 0.05

    # longwave radiation
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_lw::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_lw::FT = 0.01
end;


"""
Struct that contains the auxiliary variables of leaf biophysical traits (to save time).
"""
Base.@kwdef mutable struct LeafBioAuxil{FT<:AbstractFloat}
    # shortwave radiation
    # pigment absorption fractions
    "Chlorophyll a and b absorption fraction `[-]`"
    f_cab::Vector{FT}
    "Carotenoid absorption fraction `[-]`"
    f_car::Vector{FT}
    "PPAR fraction `[-]`"
    f_ppar::Vector{FT}
    "PSII PPAR fraction `[-]`"
    f_psii::Vector{FT}
    "SIF excitation fraction `[-]`"
    f_sife::Vector{FT}

    # reflectance and transmittance at the air-water interface
    "Air-water interface reflectance with an average angle `[-]`"
    ρ_interface_θ::Vector{FT}
    "Air-water interface transmittance with an average angle `[-]`"
    τ_interface_θ::Vector{FT}
    "Air-water interface reflectance with isotropic light `[-]`"
    ρ_interface_12::Vector{FT}
    "Air-water interface transmittance with isotropic light `[-]`"
    τ_interface_12::Vector{FT}
    "Water-air interface reflectance with isotropic light `[-]`"
    ρ_interface_21::Vector{FT}
    "Water-air interface transmittance with isotropic light `[-]`"
    τ_interface_21::Vector{FT}
    "Effective air-water interface reflectance with isotropic light `[-]`"
    ρ_interface_12_eff::Vector{FT}
    "Effective air-water interface transmittance with isotropic light `[-]`"
    τ_interface_12_eff::Vector{FT}
    "Effective water-air interface reflectance with isotropic light `[-]`"
    ρ_interface_21_eff::Vector{FT}
    "Effective water-air interface transmittance with isotropic light `[-]`"
    τ_interface_21_eff::Vector{FT}

    # transmittance within a sublayer of a layer
    "First Layer sublayer transmittance with isotropic light `[-]`"
    τ_sub_1::Vector{FT}
    "Second Layer sublayer transmittance with isotropic light `[-]`"
    τ_sub_2::Vector{FT}
    "First layer total transmittance of all sublayers `[-]`"
    τ_all_1::Vector{FT}
    "Second layer total transmittance of all sublayers `[-]`"
    τ_all_2::Vector{FT}
    "First layer total distinction coefficient of all sublayers `[-]`"
    k_all_1::Vector{FT}
    "Second layer total distinction coefficient of all sublayers `[-]`"
    k_all_2::Vector{FT}

    # reflectance and transmittance of a single layer
    "First layer reflectance with an average angle `[-]`"
    ρ_layer_θ::Vector{FT}
    "First layer transmittance with an average angle `[-]`"
    τ_layer_θ::Vector{FT}
    "First layer reflectance with isotropic light `[-]`"
    ρ_layer_1::Vector{FT}
    "First layer transmittance with isotropic light `[-]`"
    τ_layer_1::Vector{FT}
    "Second layer reflectance with isotropic light `[-]`"
    ρ_layer_2::Vector{FT}
    "Second layer transmittance with isotropic light `[-]`"
    τ_layer_2::Vector{FT}

    # reflectance and transmittance of the leaf
    "Reflectance of the leaf `[-]`"
    ρ_leaf::Vector{FT}
    "Transmittance of the leaf `[-]`"
    τ_leaf::Vector{FT}
    "Absorption of the leaf `[-]`"
    α_leaf::Vector{FT}

    # SIF excitation to emittance matrix (before scaling with Φ_PS*)
    "First layer SIF matrix backwards (emission only) `[-]`"
    mat_b_1::Matrix{FT}
    "First layer SIF matrix forwards (emission only) `[-]`"
    mat_f_1::Matrix{FT}
    "Second layer SIF matrix backwards (emission only) `[-]`"
    mat_b_2::Matrix{FT}
    "Second layer SIF matrix forwards (emission only) `[-]`"
    mat_f_2::Matrix{FT}
    "First layer SIF matrix backwards (after reabsorption, reflection, and transmission) `[-]`"
    mat_b_1_out::Matrix{FT}
    "First layer SIF matrix forwards (after reabsorption, reflection, and transmission) `[-]`"
    mat_f_1_out::Matrix{FT}
    "Second layer SIF matrix backwards (after reabsorption, reflection, and transmission) `[-]`"
    mat_b_2_out::Matrix{FT}
    "Second layer SIF matrix forwards (after reabsorption, reflection, and transmission) `[-]`"
    mat_f_2_out::Matrix{FT}

    # SIF excitation to emittance matrix
    "SIF matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "SIF matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "Mean SIF matrix of the backward and forward SIF matrices `[-]`"
    mat_mean::Matrix{FT}
    "Diff SIF matrix of the backward and forward SIF matrices `[-]`"
    mat_diff::Matrix{FT}

    # cache variables
    "SIF PDF based on the wavelength of excitation `[-]`"
    _ϕ_sif::Vector{FT}
    "SIF PDF based on the wavelength of excitation `[-]`"
    _ϕ1_sif::Vector{FT}
    "SIF PDF based on the wavelength of excitation `[-]`"
    _ϕ2_sif::Vector{FT}
end;

LeafBioAuxil(config::SPACConfig{FT}) where {FT} = (
    return LeafBioAuxil{FT}(
                f_cab              = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                f_car              = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                f_ppar             = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                f_psii             = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                f_sife             = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_interface_θ      = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_interface_θ      = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_interface_12     = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_interface_12     = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_interface_21     = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_interface_21     = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_interface_12_eff = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_interface_12_eff = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_interface_21_eff = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_interface_21_eff = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_sub_1            = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_sub_2            = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_all_1            = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_all_2            = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                k_all_1            = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                k_all_2            = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_layer_θ          = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_layer_θ          = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_layer_1          = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_layer_1          = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_layer_2          = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_layer_2          = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                ρ_leaf             = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                τ_leaf             = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                α_leaf             = zeros(FT, length(config.CONSTANTS.SPECTRA.Λ)),
                mat_b_1            = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_f_1            = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_b_2            = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_f_2            = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_b_1_out        = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_f_1_out        = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_b_2_out        = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_f_2_out        = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_b              = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_f              = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_mean           = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                mat_diff           = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF), length(config.CONSTANTS.SPECTRA.IΛ_SIFE)),
                _ϕ_sif             = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF)),
                _ϕ1_sif            = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF)),
                _ϕ2_sif            = zeros(FT, length(config.CONSTANTS.SPECTRA.IΛ_SIF)),
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: add new struct to use with the new leaf optics model (copied from HyperspectralLeafBiophysics)
#     2024-Feb-26: add field trait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafBio{FT<:AbstractFloat}
    # state variables (prognostic or structural)
    "Trait variables"
    trait::LeafBioTrait{FT} = LeafBioTrait{FT}()
    "State variables"
    state::LeafBioState{FT} = LeafBioState{FT}()
    "Auxiliary variables"
    auxil::LeafBioAuxil{FT}
end;

LeafBio(config::SPACConfig{FT}) where {FT} = return LeafBio{FT}(auxil = LeafBioAuxil(config));

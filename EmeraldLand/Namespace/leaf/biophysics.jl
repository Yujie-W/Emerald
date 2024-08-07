# This file contains the state and auxiliary variables for leaf biophysics

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct LeafBioTrait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables of leaf biophysical traits.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafBioTrait{FT<:AbstractFloat}
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

    # leaf width
    "leaf width `[m]`"
    width::FT = 0.05

    # longwave radiation
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_lw::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_lw::FT = 0.01
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: add struct LeafBioState
#     2023-Sep-14: add fields ϕ_car and ϕ_car_ppar
#     2023-Oct-03: add field width
#     2024-Jan-20: set ϕ_car_ppar to 0 by default
#     2024-Feb-26: move the traits to a separate struct
#     2024-Jan-20: set ϕ_car and ϕ_car_ppar to 1 by default
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables of leaf biophysical traits.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafBioState{FT<:AbstractFloat}
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::FT = 0
    "Fraction of carotenoid aborption into SIF `[-]`"
    ϕ_car::FT = 1
    "Fraction of carotenoid aborption into PPAR `[-]`"
    ϕ_car_ppar::FT = 1
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: add struct LeafBioAuxil
#     2023-Sep-14: add fields to store the interface, layer, and leaf reflectance and transmittance; SIF calculation matrices
#     2023-Sep-16: add fields f_sife, mat_b_i and mat_f_i
#     2023-Sep-18: add fields τ_all_i, mat_x_i_out, and _ϕ_sif
#     2023-Sep-19: add field f_ppar
#     2023-Sep-22: add field f_psii
#     2023-Oct-14: add fields mat_mean, and mat_diff
#     2023-Oct-24: add fields psi_mat_* and psii_mat_*
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables of leaf biophysical traits (to save time).

# Fields

$(TYPEDFIELDS)

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

    # transmittance within a sublayer of a layer
    "First Layer sublayer transmittance with isotropic light `[-]`"
    τ_sub_1::Vector{FT}
    "Second Layer sublayer transmittance with isotropic light `[-]`"
    τ_sub_2::Vector{FT}
    "First layer total transmittance of all sublayers `[-]`"
    τ_all_1::Vector{FT}
    "Second layer total transmittance of all sublayers `[-]`"
    τ_all_2::Vector{FT}

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

    # SIF excitation to emittance matrix for PSI
    "SIF matrix backwards `[-]`"
    psi_mat_b::Matrix{FT}
    "SIF matrix forwards `[-]`"
    psi_mat_f::Matrix{FT}
    "Mean SIF matrix of the backward and forward SIF matrices `[-]`"
    psi_mat_mean::Matrix{FT}
    "Diff SIF matrix of the backward and forward SIF matrices `[-]`"
    psi_mat_diff::Matrix{FT}

    # SIF excitation to emittance matrix for PSII
    "SIF matrix backwards `[-]`"
    psii_mat_b::Matrix{FT}
    "SIF matrix forwards `[-]`"
    psii_mat_f::Matrix{FT}
    "Mean SIF matrix of the backward and forward SIF matrices `[-]`"
    psii_mat_mean::Matrix{FT}
    "Diff SIF matrix of the backward and forward SIF matrices `[-]`"
    psii_mat_diff::Matrix{FT}

    # cache variables
    "SIF PDF based on the wavelength of excitation `[-]`"
    _ϕ_sif::Vector{FT}
    "SIF PDF based on the wavelength of excitation `[-]`"
    _ϕ1_sif::Vector{FT}
    "SIF PDF based on the wavelength of excitation `[-]`"
    _ϕ2_sif::Vector{FT}
end;

LeafBioAuxil(config::SPACConfiguration{FT}) where {FT} = (
    return LeafBioAuxil{FT}(
                f_cab            = zeros(FT, length(config.SPECTRA.Λ)),
                f_car            = zeros(FT, length(config.SPECTRA.Λ)),
                f_ppar           = zeros(FT, length(config.SPECTRA.Λ)),
                f_psii           = zeros(FT, length(config.SPECTRA.Λ)),
                f_sife           = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_interface_θ    = zeros(FT, length(config.SPECTRA.Λ)),
                τ_interface_θ    = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_interface_12   = zeros(FT, length(config.SPECTRA.Λ)),
                τ_interface_12   = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_interface_21   = zeros(FT, length(config.SPECTRA.Λ)),
                τ_interface_21   = zeros(FT, length(config.SPECTRA.Λ)),
                τ_sub_1          = zeros(FT, length(config.SPECTRA.Λ)),
                τ_sub_2          = zeros(FT, length(config.SPECTRA.Λ)),
                τ_all_1          = zeros(FT, length(config.SPECTRA.Λ)),
                τ_all_2          = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_layer_θ        = zeros(FT, length(config.SPECTRA.Λ)),
                τ_layer_θ        = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_layer_1        = zeros(FT, length(config.SPECTRA.Λ)),
                τ_layer_1        = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_layer_2        = zeros(FT, length(config.SPECTRA.Λ)),
                τ_layer_2        = zeros(FT, length(config.SPECTRA.Λ)),
                ρ_leaf           = zeros(FT, length(config.SPECTRA.Λ)),
                τ_leaf           = zeros(FT, length(config.SPECTRA.Λ)),
                α_leaf           = zeros(FT, length(config.SPECTRA.Λ)),
                mat_b_1          = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_f_1          = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_b_2          = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_f_2          = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_b_1_out      = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_f_1_out      = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_b_2_out      = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_f_2_out      = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_b            = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_f            = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_mean         = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                mat_diff         = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psi_mat_b        = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psi_mat_f        = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psi_mat_mean     = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psi_mat_diff     = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psii_mat_b       = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psii_mat_f       = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psii_mat_mean    = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                psii_mat_diff    = zeros(FT, length(config.SPECTRA.IΛ_SIF), length(config.SPECTRA.IΛ_SIFE)),
                _ϕ_sif           = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
                _ϕ1_sif          = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
                _ϕ2_sif          = zeros(FT, length(config.SPECTRA.IΛ_SIF)),
    )
);

sync_struct!(bio_from::LeafBioAuxil{FT}, bio_to::LeafBioAuxil{FT}) where {FT} = (
    bio_to.f_cab          .= bio_from.f_cab;
    bio_to.f_car          .= bio_from.f_car;
    bio_to.f_ppar         .= bio_from.f_ppar;
    bio_to.f_psii         .= bio_from.f_psii;
    bio_to.f_sife         .= bio_from.f_sife;
    bio_to.ρ_interface_θ  .= bio_from.ρ_interface_θ;
    bio_to.τ_interface_θ  .= bio_from.τ_interface_θ;
    bio_to.ρ_interface_12 .= bio_from.ρ_interface_12;
    bio_to.τ_interface_12 .= bio_from.τ_interface_12;
    bio_to.ρ_interface_21 .= bio_from.ρ_interface_21;
    bio_to.τ_interface_21 .= bio_from.τ_interface_21;
    bio_to.τ_sub_1        .= bio_from.τ_sub_1;
    bio_to.τ_sub_2        .= bio_from.τ_sub_2;
    bio_to.τ_all_1        .= bio_from.τ_all_1;
    bio_to.τ_all_2        .= bio_from.τ_all_2;
    bio_to.ρ_layer_θ      .= bio_from.ρ_layer_θ;
    bio_to.τ_layer_θ      .= bio_from.τ_layer_θ;
    bio_to.ρ_layer_1      .= bio_from.ρ_layer_1;
    bio_to.τ_layer_1      .= bio_from.τ_layer_1;
    bio_to.ρ_layer_2      .= bio_from.ρ_layer_2;
    bio_to.τ_layer_2      .= bio_from.τ_layer_2;
    bio_to.ρ_leaf         .= bio_from.ρ_leaf;
    bio_to.τ_leaf         .= bio_from.τ_leaf;
    bio_to.α_leaf         .= bio_from.α_leaf;
    bio_to.mat_b_1        .= bio_from.mat_b_1;
    bio_to.mat_f_1        .= bio_from.mat_f_1;
    bio_to.mat_b_2        .= bio_from.mat_b_2;
    bio_to.mat_f_2        .= bio_from.mat_f_2;
    bio_to.mat_b_1_out    .= bio_from.mat_b_1_out;
    bio_to.mat_f_1_out    .= bio_from.mat_f_1_out;
    bio_to.mat_b_2_out    .= bio_from.mat_b_2_out;
    bio_to.mat_f_2_out    .= bio_from.mat_f_2_out;
    bio_to.mat_b          .= bio_from.mat_b;
    bio_to.mat_f          .= bio_from.mat_f;
    bio_to.mat_mean       .= bio_from.mat_mean;
    bio_to.mat_diff       .= bio_from.mat_diff;
    bio_to.psi_mat_b      .= bio_from.psi_mat_b;
    bio_to.psi_mat_f      .= bio_from.psi_mat_f;
    bio_to.psi_mat_mean   .= bio_from.psi_mat_mean;
    bio_to.psi_mat_diff   .= bio_from.psi_mat_diff;
    bio_to.psii_mat_b     .= bio_from.psii_mat_b;
    bio_to.psii_mat_f     .= bio_from.psii_mat_f;
    bio_to.psii_mat_mean  .= bio_from.psii_mat_mean;
    bio_to.psii_mat_diff  .= bio_from.psii_mat_diff;
    bio_to._ϕ_sif         .= bio_from._ϕ_sif;
    bio_to._ϕ1_sif        .= bio_from._ϕ1_sif;
    bio_to._ϕ2_sif        .= bio_from._ϕ2_sif;

    return nothing
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

LeafBio(config::SPACConfiguration{FT}) where {FT} = return LeafBio{FT}(auxil = LeafBioAuxil(config));

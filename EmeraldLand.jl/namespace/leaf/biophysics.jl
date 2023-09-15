#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: moved the state variables out from HyperLeafBiophysics
#     2023-Sep-14: add fields ϕ_car and ϕ_car_ppar
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables of leaf biophysical traits.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperLeafBioState{FT<:AbstractFloat}
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
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::FT = 0
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
    "Leaf mesophyll structural parameter that describes the number of thin layers with a leaf"
    meso_n::FT = 1.4
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT = 0
    "Fraction of carotenoid aborption into SIF `[-]`"
    ϕ_car::FT = 0
    "Fraction of carotenoid aborption into PPAR `[-]`"
    ϕ_car_ppar::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: moved the auxiliary variables out from HyperLeafBiophysics
#     2023-Sep-14: add fields to store the interface, layer, and leaf reflectance and transmittance
#     2023-Sep-14: add fields to store the SIF calculation matrices
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables of leaf biophysical traits (to save time).

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperLeafBioAuxil{FT<:AbstractFloat}
    # shortwave radiation
    # pigment absorption fractions
    "Chlorophyll a and b absorption fraction `[-]`"
    f_cab::Vector{FT}
    "Carotenoid absorption fraction `[-]`"
    f_car::Vector{FT}

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

    # SIF excitation to emittance matrix
    "SIF matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "SIF matrix backwards without reabsorption `[-]`"
    mat_b_chl::Matrix{FT}
    "SIF matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "SIF matrix forwards without reabsorption `[-]`"
    mat_f_chl::Matrix{FT}

    # longwave radiation
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT = 0.01
end;

HyperLeafBioAuxil(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_SIF, DIM_SIFE, DIM_WL) = config;

    return HyperLeafBioAuxil{FT}(
                k_cab            = zeros(FT, DIM_WL),
                k_car            = zeros(FT, DIM_WL),
                ρ_interface_θ    = zeros(FT, DIM_WL),
                τ_interface_θ    = zeros(FT, DIM_WL),
                ρ_interface_12   = zeros(FT, DIM_WL),
                τ_interface_12   = zeros(FT, DIM_WL),
                ρ_interface_21   = zeros(FT, DIM_WL),
                τ_interface_21   = zeros(FT, DIM_WL),
                τ_sub_1          = zeros(FT, DIM_WL),
                τ_sub_2          = zeros(FT, DIM_WL),
                ρ_layer_θ        = zeros(FT, DIM_WL),
                τ_layer_θ        = zeros(FT, DIM_WL),
                ρ_layer_1        = zeros(FT, DIM_WL),
                τ_layer_1        = zeros(FT, DIM_WL),
                ρ_layer_2        = zeros(FT, DIM_WL),
                τ_layer_2        = zeros(FT, DIM_WL),
                ρ_leaf           = zeros(FT, DIM_WL),
                τ_leaf           = zeros(FT, DIM_WL),
                mat_b            = zeros(FT, DIM_SIF, DIM_SIFE),
                mat_b_chl        = zeros(FT, DIM_SIF, DIM_SIFE),
                mat_f            = zeros(FT, DIM_SIF, DIM_SIFE),
                mat_f_chl        = zeros(FT, DIM_SIF, DIM_SIFE),
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: add new struct to use with the new leaf optics model (copied from HyperspectralLeafBiophysics)
#     2023-Sep-14: clean the variables to state and auxil structs
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperLeafBio{FT<:AbstractFloat}
    # state variables (prognostic or structural)
    "State variables"
    state::HyperLeafBioState{FT} = HyperLeafBioState{FT}()
    "Auxiliary variables"
    auxil::HyperLeafBioAuxil{FT} = HyperLeafBioAuxil{FT}()
end

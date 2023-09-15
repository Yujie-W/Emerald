#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: moved the state variables out from HyperLeafBiophysics
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
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: moved the auxiliary variables out from HyperLeafBiophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables of leaf biophysical traits (to save time).

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperLeafBioAuxil{FT<:AbstractFloat}
    # longwave radiation
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT = 0.01
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-14: add new struct to use with the new leaf optics model (copied from HyperspectralLeafBiophysics)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperLeafBiophysics{FT<:AbstractFloat}
    # state variables (prognostic or structural)
    states::HyperLeafBioState{FT} = HyperLeafBioState{FT}()

    # Diagnostic variables
    "Specific absorption coefficients of all materials"
    k_all::Vector{FT}
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "Fluorescence excitation matrix backwards without reabsorption `[-]`"
    mat_b_chl::Matrix{FT}
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "Fluorescence excitation matrix forwards without reabsorption `[-]`"
    mat_f_chl::Matrix{FT}
    "Relative absorption by Chlorophyll `[-]`"
    α_cab::Vector{FT}
    "Relative absorption by Chlorophyll+Carotenoid `[-]`"
    α_cabcar::Vector{FT}
    "Shortwave absorption, 1 .- ρ_sw .- τ_sw  `[-]`"
    α_sw::Vector{FT}
    "Shortwave leaf reflectance `[-]`"
    ρ_sw::Vector{FT}
    "Shortwave leaf transmission `[-]`"
    τ_sw::Vector{FT}

    # Cache variables
    "Leaf water content history used to compute leaf spectra `[mol m⁻²]`"
    _v_storage::FT = 0

    # Cache variables to speed up leaf_spectra!
    _a::Vector{FT}
    _a²::Vector{FT}
    _b::Vector{FT}
    _bⁿ⁻¹::Vector{FT}
    _b²ⁿ⁻²::Vector{FT}
    _d::Vector{FT}
    _denom::Vector{FT}
    _k::Vector{FT}
    _k_chl::Vector{FT}
    _s::Vector{FT}
    _tt1::Vector{FT}
    _tt2::Vector{FT}
    _z::Vector{FT}

    _ρ::Vector{FT}
    _ρ²::Vector{FT}
    _ρ_b::Vector{FT}
    _ρ_bottom::Vector{FT}
    _ρ_sub::Vector{FT}
    _ρ_top::Vector{FT}
    _ρ_α::Vector{FT}
    _ρ₁₂::Vector{FT}
    _ρ₂₁::Vector{FT}
    _τ::Vector{FT}
    _τ²::Vector{FT}
    _τ_bottom::Vector{FT}
    _τ_sub::Vector{FT}
    _τ_top::Vector{FT}
    _τ_α::Vector{FT}
    _τ₁₂::Vector{FT}
    _τ₂₁::Vector{FT}

    _1_e::Matrix{FT}
    _1_f::Matrix{FT}
    _a₁₁::Matrix{FT}
    _a₁₂::Matrix{FT}
    _a₂₁::Matrix{FT}
    _a₂₂::Matrix{FT}
    _m_xe::Matrix{FT}
    _m_xf::Matrix{FT}
    _m_ye::Matrix{FT}
    _m_yf::Matrix{FT}
    _ma::Matrix{FT}
    _mb::Matrix{FT}
    _mat_b::Matrix{FT}
    _mat_b_n::Matrix{FT}
    _mat_f::Matrix{FT}
    _mat_f_n::Matrix{FT}
    _sigmoid::Matrix{FT}
    _x_e::Vector{FT}
    _x_f::Vector{FT}
    _z_e::Vector{FT}
    _z_f::Vector{FT}
    _ρ_e::Vector{FT}
    _ρ_e_n::Vector{FT}
    _ρ_f::Vector{FT}
    _ρ_f_n::Vector{FT}
    _τ_e::Vector{FT}
    _τ_e_n::Vector{FT}
    _τ_f::Vector{FT}
    _τ_f_n::Vector{FT}
end

HyperspectralLeafBiophysics(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_SIF, DIM_SIFE, DIM_WL) = config;

    return HyperspectralLeafBiophysics{FT}(
                k_all     = zeros(FT, DIM_WL),
                mat_b     = zeros(FT, DIM_SIF, DIM_SIFE),
                mat_b_chl = zeros(FT, DIM_SIF, DIM_SIFE),
                mat_f     = zeros(FT, DIM_SIF, DIM_SIFE),
                mat_f_chl = zeros(FT, DIM_SIF, DIM_SIFE),
                α_cab     = zeros(FT, DIM_WL),
                α_cabcar  = zeros(FT, DIM_WL),
                α_sw      = zeros(FT, DIM_WL),
                ρ_sw      = zeros(FT, DIM_WL),
                τ_sw      = zeros(FT, DIM_WL),
                _a        = zeros(FT, DIM_WL),
                _a²       = zeros(FT, DIM_WL),
                _b        = zeros(FT, DIM_WL),
                _bⁿ⁻¹     = zeros(FT, DIM_WL),
                _b²ⁿ⁻²    = zeros(FT, DIM_WL),
                _d        = zeros(FT, DIM_WL),
                _denom    = zeros(FT, DIM_WL),
                _k        = zeros(FT, DIM_WL),
                _k_chl    = zeros(FT, DIM_WL),
                _s        = zeros(FT, DIM_WL),
                _tt1      = zeros(FT, DIM_WL),
                _tt2      = zeros(FT, DIM_WL),
                _z        = zeros(FT, DIM_WL),
                _ρ        = zeros(FT, DIM_WL),
                _ρ²       = zeros(FT, DIM_WL),
                _ρ_b      = zeros(FT, DIM_WL),
                _ρ_bottom = zeros(FT, DIM_WL),
                _ρ_sub    = zeros(FT, DIM_WL),
                _ρ_top    = zeros(FT, DIM_WL),
                _ρ_α      = zeros(FT, DIM_WL),
                _ρ₁₂      = zeros(FT, DIM_WL),
                _ρ₂₁      = zeros(FT, DIM_WL),
                _τ        = zeros(FT, DIM_WL),
                _τ²       = zeros(FT, DIM_WL),
                _τ_bottom = zeros(FT, DIM_WL),
                _τ_sub    = zeros(FT, DIM_WL),
                _τ_top    = zeros(FT, DIM_WL),
                _τ_α      = zeros(FT, DIM_WL),
                _τ₁₂      = zeros(FT, DIM_WL),
                _τ₂₁      = zeros(FT, DIM_WL),
                _1_e      = ones(FT, 1, DIM_SIFE),
                _1_f      = ones(FT, DIM_SIF, 1),
                _a₁₁      = zeros(FT, DIM_SIF, DIM_SIFE),
                _a₁₂      = zeros(FT, DIM_SIF, DIM_SIFE),
                _a₂₁      = zeros(FT, DIM_SIF, DIM_SIFE),
                _a₂₂      = zeros(FT, DIM_SIF, DIM_SIFE),
                _m_xe     = zeros(FT, DIM_SIF, DIM_SIFE),
                _m_xf     = zeros(FT, DIM_SIF, DIM_SIFE),
                _m_ye     = zeros(FT, DIM_SIF, DIM_SIFE),
                _m_yf     = zeros(FT, DIM_SIF, DIM_SIFE),
                _ma       = zeros(FT, DIM_SIF, DIM_SIFE),
                _mb       = zeros(FT, DIM_SIF, DIM_SIFE),
                _mat_b    = zeros(FT, DIM_SIF, DIM_SIFE),
                _mat_b_n  = zeros(FT, DIM_SIF, DIM_SIFE),
                _mat_f    = zeros(FT, DIM_SIF, DIM_SIFE),
                _mat_f_n  = zeros(FT, DIM_SIF, DIM_SIFE),
                _sigmoid  = zeros(FT, DIM_SIF, DIM_SIFE),
                _x_e      = zeros(FT, DIM_SIFE),
                _x_f      = zeros(FT, DIM_SIF),
                _z_e      = zeros(FT, DIM_SIFE),
                _z_f      = zeros(FT, DIM_SIF),
                _ρ_e      = zeros(FT, DIM_SIFE),
                _ρ_e_n    = zeros(FT, DIM_SIFE),
                _ρ_f      = zeros(FT, DIM_SIF),
                _ρ_f_n    = zeros(FT, DIM_SIF),
                _τ_e      = zeros(FT, DIM_SIFE),
                _τ_e_n    = zeros(FT, DIM_SIFE),
                _τ_f      = zeros(FT, DIM_SIF),
                _τ_f_n    = zeros(FT, DIM_SIF),
    )
);

#=

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for leaf biophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractLeafBiophysics:
- [`BroadbandLeafBiophysics`](@ref)
- [`HyperspectralLeafBiophysics`](@ref)

"""
abstract type AbstractLeafBiophysics{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-15: add struct for broadband leaf biophysics
#     2022-Jun-24: add leaf emissivity constant
#     2022-Jul-19: add field lma
#     2022-Jul-20: use α and ϵ instead of upper case greek Α and Ε
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandLeafBiophysics{FT<:AbstractFloat} <: AbstractLeafBiophysics{FT}
    # General information of leaf biophysics
    "Broadband absorption fraction at the NIR region"
    α_NIR::FT = 0.2
    "Broadband absorption fraction at the PAR region"
    α_PAR::FT = 0.8
    "Emissivity for longwave radiation"
    ϵ_LW::FT = 0.97

    # Prognostic variables
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jun-15: rename struct to HyperspectralLeafBiophysics to distinguish from BroadbandLeafBiophysics
#     2022-Jul-19: add dimension control to struct
#     2022-Jul-20: remove field: l_H₂O (unit cm)
#     2022-Jul-20: rename ρ_lw and τ_lw to ρ_LW and τ_LW
#     2022-Jul-28: add field _v_storage to speed up calculations (run leaf_spectra! only of _v_storage differs from current leaf water content)
#     2023-Jun-16: remove fields DIM_*
#     2023-Sep-09: add field mat_b_chl and mat_f_chl
#     2023-Sep-11: add mat_b_chl and mat_f_chl in the constructor function
#     2023-Sep-09: remove field mat_b_chl and mat_f_chl because SCOPE scheme is incorrect
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralLeafBiophysics{FT<:AbstractFloat} <: AbstractLeafBiophysics{FT}
    # General information of leaf biophysics
    "Leaf mesophyll structural parameter that describes mesophyll reflectance and transmittance"
    MESOPHYLL_N::FT = 1.4
    "Doubling adding layers"
    NDUB::Int = 10
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT = 0.01

    # Prognostic variables
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
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT = 0

    # Diagnostic variables
    "Specific absorption coefficients of all materials"
    k_all::Vector{FT}
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT}
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
                mat_f     = zeros(FT, DIM_SIF, DIM_SIFE),
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

=#


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2021-Nov-29: add abstract photosynthesis system type
#     2022-Jan-14: rename to photosynthesis model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractPhotosynthesisModel`:
- [`C3CytochromeModel`](@ref)
- [`C3VJPModel`](@ref)
- [`C4VJPModel`](@ref)

"""
abstract type AbstractPhotosynthesisModel{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add C3CytochromeModel structure for C₃ photosynthesis system
#     2022-Feb-07: add more fields to use with Photosynthesis v0.3.1
#     2022-Feb-07: remove j_p680 and j_p700 series variables
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Mar-01: add two more fields: TD_ΗC and TD_ΗL
#     2022-Mar-01: move η_c and η_l from reaction center to photosynthesis model
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: remove v_cmax25_ww (use β instead)
#     2022-Jul-20: rename TD_ΗC and TD_ΗL to lower greek TD_ηC and TD_ηL
#     2022-Jul-20: move a_c, a_j, and etc as cache variables, and rename them to _a_c, _a_j, and etc
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C3 Cytochrome photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3CytochromeModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = SerialColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KcTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KoTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kq temperature dependency"
    TD_KQ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KqTDJohnson(FT)
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ΓStarTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_C temperature dependency"
    TD_ηC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηCTDJohnson(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_L temperature dependency"
    TD_ηL::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηLTDJohnson(FT)

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_c::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_j::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_p::FT = 0
    "Electron to CO₂ coefficient"
    _e_to_c::FT = 0
    "Potential Electron Transport Rate `[μmol e⁻ m⁻² s⁻¹]`"
    _j_pot::FT = 0
    "PSI electron transport rate after colimitation"
    _j_psi::FT = 0
    "RubisCO coefficient Kc `[Pa]`"
    _k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    _k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    _k_o::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    _k_q::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _r_d::FT = 0
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_cmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    _v_qmax::FT = 0
    "ratio between J_P700 and J_P680"
    _η::FT = 0
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    _η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    _η_l::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    _γ_star::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Nov-11: add C3VJPModel structure for classic C₃ photosynthesis system
#     2022-Jan-14: add temperature dependency into the structure
#     2022-Jan-14: rename to photosynthesis model
#     2022-Jan-14: add colimitation and e_to_c
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ, COLIMIT_IP, and COLIMIT_J (minor breaking)
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: remove v_cmax25_ww (use β instead)
#     2022-Jul-20: move a_c, a_j, and etc as cache variables, and rename them to _a_c, _a_j, and etc
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C3 photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3VJPModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = QuadraticColimit{FT}(CURVATURE = 0.7)

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Jmax temperature dependency"
    TD_JMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = JmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KcTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KoTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ΓStarTDCLM(FT)

    # Prognostic variables
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT = 83.5
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_c::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_j::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_p::FT = 0
    "Electron to CO₂ coefficient"
    _e_to_c::FT = 0
    "Electron transport `[μmol m⁻² s⁻¹]`"
    _j::FT = 0
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _j_max::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    _j_pot::FT = 0
    "RubisCO coefficient Kc `[Pa]`"
    _k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    _k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    _k_o::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _r_d::FT = 0
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_cmax::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    _γ_star::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add C4VJPModel structure for classic C₄ photosynthesis system
#     2022-Feb-11: remove j from the struct
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: remove v_cmax25_ww (use β instead)
#     2022-Jul-20: move a_c, a_j, and etc as cache variables, and rename them to _a_c, _a_j, and etc
# To do
#     TODO: add Jmax to C4VJPModel and thus JMAX TD in Photosynthesis.jl (not necessary)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C4 photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4VJPModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency"
    TD_KPEP::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KpepTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type  respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vpmax temperature dependency"
    TD_VPMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VpmaxTDBoyd(FT)

    # Prognostic variables
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT = 50

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_c::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_j::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_p::FT = 0
    "Electron to CO₂ coefficient"
    _e_to_c::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    _j_pot::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    _k_pep::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _r_d::FT = 0
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_pmax::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-14: add abstract mode type
#     2022-Jan-14: add Hierarchy description
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractPhotosynthesisMode:
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)

"""
abstract type AbstractPhotosynthesisMode end


""" An empty structure to signal the function to calculate photosynthetic rates based on leaf diffusive conductance to CO₂ """
struct GCO₂Mode <: AbstractPhotosynthesisMode end


""" An empty structure to signal the function to calculate photosynthetic rates based on CO₂ partial pressure """
struct PCO₂Mode <: AbstractPhotosynthesisMode end


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2022-Jan-14: add van der Tol model struct
#     2022-Feb-07: remove the Hierarchy from abstract fluorescence model
#     2022-Feb-07: move struct definition as a field of VJPReactionCenter
#     2022-Oct-27: make VanDerTolFluorescenceModel struct mutable
# To do
#     TODO: examine why van der Tol et al has the nonstressed parameter set that are so far away from the stressed one
# Sources
#     van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores van der Tol et al. (2014) fluorescence model parameters.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VanDerTolFluorescenceModel{FT<:AbstractFloat}
    # General model information
    "Fitting parameter K_0"
    K_0::FT = 5.01
    "Fitting parameter α"
    K_A::FT = 1.93
    "Fitting parameter β"
    K_B::FT = 10
end


""" VanDerTolFluorescenceModel that uses data from all observations """
VDTModelAll(FT) = VanDerTolFluorescenceModel{FT}(K_0 = 2.48, K_A = 2.83, K_B = 0.114)


""" VanDerTolFluorescenceModel that uses data from drought stressed observations """
VDTModelDrought(FT) = VanDerTolFluorescenceModel{FT}(K_0 = 5.01, K_A = 1.93, K_B = 10);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-24: abstractize the reaction center
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Abstract type for reaction center

Hierarchy of the `AbstractReactionCenter`
- [`VJPReactionCenter`](@ref)
- [`CytochromeReactionCenter`](@ref)

"""
abstract type AbstractReactionCenter{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-15: isolate the reaction center from Leaf in Photosynthesis.jl
#     2022-Feb-07: add fluorescence model as a field
#     2022-Jul-20: rename quite some cache variables
#     2023-Jan-27: set default K_D to 0.95 instead of 0.85
#     2023-Sep-09: add fields ϕ_d and ϕ_n
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VJPReactionCenter{FT<:AbstractFloat} <:AbstractReactionCenter{FT}
    # Constant coefficients
    "Fraction of absorbed light used by PSII ETR"
    F_PSII::FT = 0.5
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for photochemistry"
    K_P_MAX::FT = 4

    # Embedded structures
    "Fluorescence model"
    FLM::VanDerTolFluorescenceModel{FT} = VanDerTolFluorescenceModel{FT}()

    # Prognostic variables
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0

    # Diagnostic variables
    "Heat dissipation yield"
    ϕ_d::FT = 0
    "Fluorescence yield"
    ϕ_f::FT = 0
    "Non-photochemical quenching yeild"
    ϕ_n::FT = 0
    "Photochemical yield"
    ϕ_p::FT = 0

    # Cache variables
    "Dark adapted yield (`Kp=0`)"
    _f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    _f_m′::FT = 0
    "Dark-adapted fluorescence yield (`Kp=max`)"
    _f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    _f_o′::FT = 0
    "Rate constant for thermal dissipation"
    _k_d::FT = 0.95
    "Reversible NPQ rate constant (initially zero)"
    _k_npq_rev::FT = 0
    "Rate constant for photochemistry"
    _k_p::FT = 4
    "Non-Photochemical quenching "
    _npq::FT = 0
    "Energy quenching"
    _q_e::FT = 0
    "Photochemical quenching"
    _q_p::FT = 0
    "max PSII yield (_k_npq_rev = 0, all RC open)"
    _ϕ_psii_max::FT = K_P_MAX / (_k_d + K_F + K_P_MAX)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add the struct of Cytochrome reaction center
#     2022-Feb-07: add more fields
#     2022-Feb-10: add K_X, ϵ_1, and ϵ_2 fields
#     2022-Mar-01: add more fields: η_c and η_l
#     2022-Mar-01: delete Η_C and Η_L, move η_c and η_l to photosynthesis model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CytochromeReactionCenter{FT<:AbstractFloat} <:AbstractReactionCenter{FT}
    # Constant coefficients
    "Fraction of absorbed light used by PSI ETR"
    F_PSI::FT = 0.41 / (0.41 + 0.44)
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.55
    "Rate constant of fluorescence `[ns⁻¹]`"
    K_F::FT = 0.05
    "Rate constant of photochemistry for PS I `[ns⁻¹]`"
    K_PSI::FT = 14.5
    "Rate constant of photochemistry for PS II `[ns⁻¹]`"
    K_PSII::FT = 4.5
    "Rate constant of excitation sharing for PS II `[ns⁻¹]`"
    K_U::FT = 2
    "Rate constant of regulated heat loss via oxidized PS I center `[s⁻¹]`"
    K_X::FT = 14.5
    "Maximal PS I photochemical yield"
    Φ_PSI_MAX::FT = K_PSI / (K_D + K_F + K_PSI)

    # Diagnostic variables
    "Weight factor that PSI fluorescence reaches sensor (after reabsorption)"
    ϵ_1::FT = 0
    "Weight factor that PSII fluorescence reaches sensor (after reabsorption)"
    ϵ_2::FT = 1
    "Fluorescence yield"
    ϕ_f::FT = 0
    "Photochemical yield"
    ϕ_p::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jul-19: abstractize the leaf
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Abstract type for leaf

Hierarchy of the `AbstractLeaf`
- [`Leaf`](@ref)
- [`Leaves1D`](@ref)
- [`Leaves2D`](@ref)

"""
abstract type AbstractLeaf{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jan-14: refactor the Leaf structure within BIO, PRC, PSM as fields
#     2022-Jan-24: add p_CO₂_s to the structure
#     2022-Feb-07: moved FLM to PRC
#     2022-May-25: add new field HS, WIDTH
#     2022-Jun-14: use Union instead of Abstract... for type definition
#     2022-Jun-15: add support to BroadbandLeafBiophysics and HyperspectralLeafBiophysics types
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add fields: G_LIMITS, a_gross and a_net
#     2022-Jul-12: add field: ∂g∂t
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Jul-28: move field _t to PSM
#     2022-Nov-18: use Union type for SM
#     2023-Mar-02: set minimum G to 1e-4 instead of 1e-2
#     2023-Apr-13: move field APAR_CAR to SPACConfiguration
#     2023-Jun-13: add field: etr
#     2023-Sep-11: set minimum G to 0 instead of 1e-4
#     2023-Sep-18: use HyperLeafBio instead of HyperspectralLeafBiophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters. This structure is meant for leaf level research and canopy radiative transfer scheme without sunlit and shaded partitioning (ppar and ppar-dependent variables).

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaf2{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # Constants
    # "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    # CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-3, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "New leaf struct, will replace Leaf2 in the next major refactor"
    NS::Leaf{FT} = Leaf{FT}()
    # "[`AbstractLeafBiophysics`](@ref) type leaf biophysical parameters"
    # BIO::Union{BroadbandLeafBiophysics{FT}, HyperLeafBio{FT}}
    # "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    # HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::FT = 1000
    # "Current leaf temperature"
    # t::FT = T₂₅(FT)

    # Prognostic variables (used for ∂y∂t)
    # "Total stored energy per area `[J m⁻²]`"
    # e::FT = (CP * BIO.state.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * t
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::FT = 0.01
    # "Marginal increase in energy `[W m⁻²]`"
    # ∂e∂t::FT = 0
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::FT = 0

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0
    "Actual electron transport `[μmol m⁻² s⁻¹]`"
    etr::FT = 0

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::FT = 0
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    _g_CO₂::FT = 0
    "Leaf internal CO₂ partial pressure `[Pa]`"
    _p_CO₂_i::FT = 0
    "Leaf surface CO₂ partial pressure `[Pa]`"
    _p_CO₂_s::FT = 0
end

Leaf2(config::SPACConfiguration{FT}) where {FT} = (
    return Leaf2{FT}(
                NS = Leaf(config)
    )
)


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-27: add new structure for leaves with 1D Vector of parameters, such as leaves for sunlit and shaded partitions
#     2022-Jun-27: make BIO BroadbandLeafBiophysics only
#     2022-Jun-28: add a_gross and a_net, make t a Vector, remove _t
#     2022-Jun-30: add a second HS2 for shaded leaves
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#     2022-Jul-07: make p_H₂O_sat a vector
#     2022-Jul-12: add field: ∂g∂t
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Nov-18: use Union type for SM
#     2023-Mar-02: set minimum G to 1e-4 instead of 1e-2
#     2023-Jun-13: add field: etr
#     2023-Sep-11: set minimum G to 0 instead of 1e-4
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters for a single canopy layer. This structure is meant for canopy level research and canopy radiative transfer scheme with sunlit and shaded partitioning.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaves1D{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # Constants
    # "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    # CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-3, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "New leaf struct, will replace Leaf2 in the next major refactor"
    NS::Leaf{FT}
    # "[`BroadbandLeafBiophysics`](@ref) type leaf biophysical parameters"
    # BIO::BroadbandLeafBiophysics{FT} = BroadbandLeafBiophysics{FT}()
    # "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    # HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    # "[`LeafHydraulics`](@ref) type leaf hydraulic system used for other calculations (say sunlit and shaded leaf partitioning)"
    # HS2::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::Vector{FT} = FT[3, 3]
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::Vector{FT} = FT[1000, 200]
    # "Current leaf temperature"
    # t::Vector{FT} = FT[T₂₅(FT), T₂₅(FT)]

    # Prognostic variables (used for ∂y∂t)
    # "Total stored energy per area `[J m⁻²]`"
    # e::Vector{FT} = FT[(CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * t[1], (CP * BIO.lma * 10 + HS2.v_storage * CP_L_MOL(FT)) * t[2]]
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::Vector{FT} = FT[0.01, 0.01]
    # "Marginal increase in energy `[W m⁻²]`"
    # ∂e∂t::Vector{FT} = FT[0, 0]
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::Vector{FT} = FT[0, 0]

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::Vector{FT} = FT[0, 0]
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::Vector{FT} = FT[0, 0]
    "Actual electron transport `[μmol m⁻² s⁻¹]`"
    etr::Vector{FT} = FT[0, 0]

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::Vector{FT} = FT[0, 0]
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    _g_CO₂::Vector{FT} = FT[0, 0]
    "Leaf internal CO₂ partial pressure `[Pa]`"
    _p_CO₂_i::Vector{FT} = FT[0, 0]
    "Leaf surface CO₂ partial pressure `[Pa]`"
    _p_CO₂_s::Vector{FT} = FT[0, 0]
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-27: add new structure for leaves with 2D Matrix of parameters for sunlit partitioning and point value for shaded partitioning
#     2022-Jun-27: make BIO HyperspectralLeafBiophysics only
#     2022-Jun-27: add sunlit and shaded ppar to struct (remove the ppar in canopy radiation)
#     2022-Jun-28: add a_gross, a_net, and ϕ_f for sunlit and shaded leaves
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#     2022-Jul-12: add fields: ∂g∂t_shaded and ∂g∂t_sunlit
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Jul-19: add dimension control to struct
#     2022-Jul-28: move field _t to PSM
#     2022-Nov-18: use Union type for SM
#     2023-Mar-02: set minimum G to 1e-4 instead of 1e-2
#     2023-Apr-13: move field APAR_CAR to SPACConfiguration
#     2023-Jun-13: add field: etr_shaded, etr_sunlit
#     2023-Jun-16: remove fields DIM_*
#     2023-Sep-07: add water flow integrators
#     2023-Sep-09: add fields ϕ_x_shaded and ϕ_x_sunlit
#     2023-Sep-11: set minimum G to 0 instead of 1e-4
#     2023-Sep-18: use HyperLeafBio instead of HyperspectralLeafBiophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters for a single canopy layer. This structure is meant for canopy level research and canopy radiative transfer scheme with sunlit and shaded partitioning as well as leaf
    angular distribution.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaves2D{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # Constants
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-3, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "New leaf struct, will replace Leaf2 in the next major refactor"
    NS::Leaf{FT}
    # "[`HyperLeafBio`](@ref) type leaf biophysical parameters"
    # BIO::HyperLeafBio{FT}
    # "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    # HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Absorbed photosynthetically active radiation used for photosynthesis for shaded leaves `[μmol m⁻² s⁻¹]`"
    ppar_shaded::FT = 200
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit leaves `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Matrix{FT} =
    # "Current leaf temperature `[K]`"
    # t::FT = T₂₅(FT)

    # Prognostic variables (used for ∂y∂t)
    # "Total stored energy per area `[J m⁻²]`"
    # e::FT = (CP * BIO.state.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * t
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT = 0.01
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT}
    # "Marginal increase in energy `[W m⁻²]`"
    # ∂e∂t::FT = 0
    "Marginal increase of conductance per time for shaded leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_shaded::FT = 0
    "Marginal increase of conductance per timefor sunlit leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_sunlit::Matrix{FT}

    # Diagnostic variables
    "Gross photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_gross_shaded::FT = 0
    "Gross photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_gross_sunlit::Matrix{FT}
    "Net photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_net_shaded::FT = 0
    "Net photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_net_sunlit::Matrix{FT}
    "Actual electron transport for shaded leaves `[μmol m⁻² s⁻¹]`"
    etr_shaded::FT = 0
    "Actual electron transport for sunlit leaves `[μmol m⁻² s⁻¹]`"
    etr_sunlit::Matrix{FT}
    "Heat dissipation quantum yield for shaded leaves `[-]`"
    ϕ_d_shaded::FT = 0
    "Heat dissipation quantum yield for sunlit leaves `[-]`"
    ϕ_d_sunlit::Matrix{FT}
    "Fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f_shaded::FT = 0
    "Fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f_sunlit::Matrix{FT}
    "Non-photochemical quenching quantum yield for shaded leaves `[-]`"
    ϕ_n_shaded::FT = 0
    "Non-photochemical quenching quantum yield for sunlit leaves `[-]`"
    ϕ_n_sunlit::Matrix{FT}
    "Photochemical quantum yield for shaded leaves `[-]`"
    ϕ_p_shaded::FT = 0
    "Photochemical quantum yield for sunlit leaves `[-]`"
    ϕ_p_sunlit::Matrix{FT}
    "Integrator for transpiration in"
    ∫∂w∂t_in = 0
    "Integrator for transpiration out"
    ∫∂w∂t_out = 0

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::FT = 0
    "Total leaf diffusive conductance to CO₂ for shaded leaves `[mol m⁻² s⁻¹]`"
    _g_CO₂_shaded::FT = 0
    "Total leaf diffusive conductance to CO₂ for sunlit leaves `[mol m⁻² s⁻¹]`"
    _g_CO₂_sunlit::Matrix{FT}
    "Leaf internal CO₂ partial pressure for shaded leaves `[Pa]`"
    _p_CO₂_i_shaded::FT = 0
    "Leaf internal CO₂ partial pressure for sunlit leaves `[Pa]`"
    _p_CO₂_i_sunlit::Matrix{FT}
    "Leaf surface CO₂ partial pressure for shaded leaves `[Pa]`"
    _p_CO₂_s_shaded::FT = 0
    "Leaf surface CO₂ partial pressure for sunlit leaves `[Pa]`"
    _p_CO₂_s_sunlit::Matrix{FT}
end

Leaves2D(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_AZI, DIM_INCL) = config;

    return Leaves2D{FT}(
                BIO             = HyperLeafBio(config),
                ppar_sunlit     = 1000 .* ones(FT, DIM_INCL, DIM_AZI),
                g_H₂O_s_sunlit  = FT(0.01) .* ones(FT, DIM_INCL, DIM_AZI),
                ∂g∂t_sunlit     = zeros(FT, DIM_INCL, DIM_AZI),
                a_gross_sunlit  = zeros(FT, DIM_INCL, DIM_AZI),
                a_net_sunlit    = zeros(FT, DIM_INCL, DIM_AZI),
                etr_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_d_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_f_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_n_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                ϕ_p_sunlit      = zeros(FT, DIM_INCL, DIM_AZI),
                _g_CO₂_sunlit   = zeros(FT, DIM_INCL, DIM_AZI),
                _p_CO₂_i_sunlit = zeros(FT, DIM_INCL, DIM_AZI),
                _p_CO₂_s_sunlit = zeros(FT, DIM_INCL, DIM_AZI),
    )
);

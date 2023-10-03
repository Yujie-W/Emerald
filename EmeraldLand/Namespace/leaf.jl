

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
Base.@kwdef mutable struct Leaves2D{FT}
    # Constants
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-3, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "New leaf struct, will replace Leaves2D in the next major refactor"
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

    # Prognostic variables (used for ∂y∂t)
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT = 0.01
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT}
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
                NS              = Leaf(config),
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

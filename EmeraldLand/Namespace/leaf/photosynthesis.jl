# This file contains the state and auxiliary variables for leaf photosynthesis

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-27: use modified TD for η_C and η_L
#     2024-Jul-30: add K_OCS to compute internal conductance for OCS
#     2024-Jul-31: add new GeneralC3Trait struct
#     2024-Aug-01: add support to Q10Peak, Q10PeakHT, and Q10PeakLTHT
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GeneralC3Trait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = ColimitJCLM(FT)

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Jmax temperature dependency"
    TD_JMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = JmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = KcTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = KoTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kq temperature dependency"
    TD_KQ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = KqTDJohnson(FT)
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = RespirationTDCLMC3(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = VcmaxTDCLMC3(FT)
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = ΓStarTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_C temperature dependency"
    TD_ηC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = ηCTDWang(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_L temperature dependency"
    TD_ηL::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = ηLTDWang(FT)

    # Related to OCS uptake
    "Multiplier to derive internal conductance for OCS `[mol μmol⁻¹]`"
    K_OCS::FT = 1400 * 1e-6

    # Embedded method structures
    "Ac method"
    ACM::AcMethodC3VcmaxPi = AcMethodC3VcmaxPi()
    "Aj method"
    AJM::Union{AjMethodC3JmaxPi, AjMethodC3VqmaxPi} = AjMethodC3JmaxPi()
    "Ap method"
    APM::Union{ApMethodC3Inf, ApMethodC3Vcmax} = ApMethodC3Vcmax()
    "Fluorescence model"
    FLM::Union{CytochromeFluoscenceModel{FT}, KNFluoscenceModel{FT}, QLFluoscenceModel{FT}, QLFluoscenceModelHan{FT}} = KNFluoscenceModel{FT}()

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT = 83.5
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Apr-15: add C4CLMTrait struct
#     2024-Jul-30: add K_OCS to compute internal conductance for OCS
#     2024-Jul-31: add new GeneralC4Trait struct
#     2024-Aug-01: add support to Q10Peak, Q10PeakHT, and Q10PeakLTHT
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C4 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GeneralC4Trait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency"
    TD_KPEP::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = KpepTDBoyd(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency to use with C4CLM method"
    TD_KPEP_CLM::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = Q10TDKpepCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type  respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = RespirationTDCLMC4(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = VcmaxTDCLMC4(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vpmax temperature dependency"
    TD_VPMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}} = VpmaxTDBoyd(FT)

    # Related to OCS uptake
    "Multiplier to derive internal conductance for OCS `[mol μmol⁻¹]`"
    K_OCS::FT = 8862 * 1e-6

    # Embedded structures
    "Ac method"
    ACM::AcMethodC4Vcmax = AcMethodC4Vcmax()
    "Aj method"
    AJM::AjMethodC4JPSII = AjMethodC4JPSII()
    "Ap method"
    APM::Union{ApMethodC4VcmaxPi, ApMethodC4VpmaxPi} = ApMethodC4VcmaxPi()
    "Fluorescence model"
    FLM::Union{KNFluoscenceModel{FT}, QLFluoscenceModel{FT}, QLFluoscenceModelHan{FT}} = KNFluoscenceModel{FT}()

    # Prognostic variables
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C3State struct
#     2023-Oct-28: add support to QLFluoscenceModel
#     2024-Jul-22: support all C3 models
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3State{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Prognostic variables (for VJP model)
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C4State struct
#     2023-Oct-28: add support to QLFluoscenceModel
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for C4 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4State{FT}
    # Prognostic variables
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add LeafPhotosystemAuxil struct
#     2023-Oct-24: add fields ϕ_f1 and ϕ_f2; remove fields ϵ_1 and ϵ_2 (computed in the LeafOptics module)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables for leaf photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafPhotosystemAuxil{FT}
    # photosynthetic rates
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_c::FT = 0
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_g::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_j::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_n::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_p::FT = 0

    # electron transport rates
    "Electron to CO₂ coefficient"
    e2c::FT = 0
    "Fraction of absorbed light used by PSII ETR"
    f_psii::FT = 0.5
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::FT = 0
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    j_max::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    j_pot::FT = 0
    "PSI electron transport rate after colimitation"
    j_psi::FT = 0

    # photosynthesis rate coefficients
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    k_pep::FT = 0
    "PEP coefficient Kpep fro CLM (different algorithm) `[Pa]`"
    k_pep_clm::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT = 0

    # respiration and carboxylation
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_pmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT = 0

    # C3 Cytochrome model variables
    "ratio between J_P700 and J_P680"
    η::FT = 0
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT = 0

    # yield variables
    "Fluorescence yield"
    ϕ_f::FT = 0
    "Photochemical yield"
    ϕ_p::FT = 0

    # fluorescence yeolds of two photosystems
    "Fluorescence yield of PSI"
    ϕ_f1::FT = 0
    "Fluorescence yield of PSII"
    ϕ_f2::FT = 0

    # fluorescence variables
    "Dark adapted yield (`Kp=0`)"
    f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    f_m′::FT = 0
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    f_o′::FT = 0
    "Non-Photochemical quenching "
    npq::FT = 0
    "Energy quenching"
    q_e::FT = 0
    "Photochemical quenching"
    q_p::FT = 0

    # fluorescence rate coefficients
    "Rate constant for thermal dissipation"
    k_d::FT = 0
    "Reversible NPQ rate constant (initially zero)"
    k_n::FT = 0
    "Rate constant for photochemistry"
    k_p::FT = 0
    "Maximal PS I photochemical yield"
    ϕ_psi_max::FT = 0
    "max PSII yield (_k_npq_rev = 0, all RC open)"
    ϕ_psii_max::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-25: define CanopyLayerPhotosystemAuxil struct to store 1D leaf photosynthesis variables (for canopy layer; Leaf will be repurposed back to elementwise)
#     2024-Jul-30: do not bin PPAR if DIM_PPAR_BINS is nothing
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables for leaf photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyLayerPhotosystemAuxil{FT}
    # photosynthetic rates
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_c::Vector{FT}
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_g::Vector{FT}
    "Intermediate photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_i::Vector{FT}
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_j::Vector{FT}
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_n::Vector{FT}
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_p::Vector{FT}

    # electron transport rates
    "Electron to CO₂ coefficient"
    e2c::Vector{FT}
    "Fraction of absorbed light used by PSII ETR"
    f_psii::FT = 0.5
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::Vector{FT}
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    j_max::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    j_pot::Vector{FT}
    "PSI electron transport rate after colimitation"
    j_psi::Vector{FT}

    # photosynthesis rate coefficients
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    k_pep::FT = 0
    "PEP coefficient Kpep fro CLM (different algorithm) `[Pa]`"
    k_pep_clm::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT = 0

    # respiration and carboxylation
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_pmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT = 0

    # C3 Cytochrome model variables
    "ratio between J_P700 and J_P680"
    η::Vector{FT}
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT = 0

    # yield variables
    "Fluorescence yield"
    ϕ_f::Vector{FT}
    "Photochemical yield"
    ϕ_p::Vector{FT}

    # fluorescence yeolds of two photosystems
    "Fluorescence yield of PSI"
    ϕ_f1::Vector{FT}
    "Fluorescence yield of PSII"
    ϕ_f2::Vector{FT}

    # fluorescence variables
    "Dark adapted yield (`Kp=0`)"
    f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    f_m′::Vector{FT}
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    f_o′::Vector{FT}
    "Non-Photochemical quenching "
    npq::Vector{FT}
    "Energy quenching"
    q_e::Vector{FT}
    "Photochemical quenching"
    q_p::Vector{FT}

    # fluorescence rate coefficients
    "Rate constant for thermal dissipation"
    k_d::FT = 0
    "Reversible NPQ rate constant (initially zero)"
    k_n::Vector{FT}
    "Rate constant for photochemistry"
    k_p::Vector{FT}
    "Maximal PS I photochemical yield"
    ϕ_psi_max::FT = 0
    "max PSII yield (_k_npq_rev = 0, all RC open)"
    ϕ_psii_max::FT = 0
end;

CanopyLayerPhotosystemAuxil(config::SPACConfiguration{FT}) where {FT} = (
    cache_dim_ppar = isnothing(config.DIM_PPAR_BINS) ? config.DIM_INCL * config.DIM_AZI : config.DIM_PPAR_BINS;

    return CanopyLayerPhotosystemAuxil{FT}(
                a_c   = zeros(FT, cache_dim_ppar+1),
                a_g   = zeros(FT, cache_dim_ppar+1),
                a_i   = zeros(FT, cache_dim_ppar+1),
                a_j   = zeros(FT, cache_dim_ppar+1),
                a_n   = zeros(FT, cache_dim_ppar+1),
                a_p   = zeros(FT, cache_dim_ppar+1),
                e2c   = zeros(FT, cache_dim_ppar+1),
                j     = zeros(FT, cache_dim_ppar+1),
                j_pot = zeros(FT, cache_dim_ppar+1),
                j_psi = zeros(FT, cache_dim_ppar+1),
                η     = zeros(FT, cache_dim_ppar+1),
                ϕ_f   = zeros(FT, cache_dim_ppar+1),
                ϕ_p   = zeros(FT, cache_dim_ppar+1),
                ϕ_f1  = zeros(FT, cache_dim_ppar+1),
                ϕ_f2  = zeros(FT, cache_dim_ppar+1),
                f_m′  = zeros(FT, cache_dim_ppar+1),
                f_o′  = zeros(FT, cache_dim_ppar+1),
                npq   = zeros(FT, cache_dim_ppar+1),
                q_e   = zeros(FT, cache_dim_ppar+1),
                q_p   = zeros(FT, cache_dim_ppar+1),
                k_n   = zeros(FT, cache_dim_ppar+1),
                k_p   = zeros(FT, cache_dim_ppar+1)
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C3VJP, C3Cyto, and C4VJP structs
#     2023-Oct-36: combine C3Cyto, C3VJP, and C4VJP into LeafPhotosystem
#     2024-Feb-26: add field trait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the fields for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafPhotosystem{FT}
    "Trait variables"
    trait::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}} = GeneralC3Trait{FT}()
    "State variables"
    state::Union{C3State{FT}, C4State{FT}} = C3State{FT}()
    "Auxilary variables"
    auxil::LeafPhotosystemAuxil{FT} = LeafPhotosystemAuxil{FT}()
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-25: add CanopyLayerPhotosystem
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the fields for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CanopyLayerPhotosystem{FT}
    "Trait variables"
    trait::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}} = GeneralC3Trait{FT}()
    "State variables"
    state::Union{C3State{FT}, C4State{FT}} = C3State{FT}()
    "Auxilary variables"
    auxil::CanopyLayerPhotosystemAuxil{FT}
end;

CanopyLayerPhotosystem(config::SPACConfiguration{FT}) where {FT} = CanopyLayerPhotosystem{FT}(auxil = CanopyLayerPhotosystemAuxil(config));

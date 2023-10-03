# This file contains the state and auxilary variables for leaf photosynthesis

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C3VJPState struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3VJPState{FT}
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
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT = 83.5
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C3CytoState struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for C3 photosynthesis (Cytochrome model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3CytoState{FT}
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
    TD_ηL::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηLTDJohnson(FT)\

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

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C4VJPState struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for C4 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4VJPState{FT}
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
#     2023-Oct-03: add PSMAuxil struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxilary variables for leaf photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PSMAuxil{FT}
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
    "Electron to CO₂ coefficient"
    e2c::FT = 0
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::FT = 0
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    j_max::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    j_pot::FT = 0
    "PSI electron transport rate after colimitation"
    j_psi::FT = 0
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    k_pep::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_pmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT = 0
    "ratio between J_P700 and J_P680"
    η::FT = 0
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT = 0

    # VJP reaction center vars
    "Heat dissipation yield"
    ϕ_d::FT = 0
    "Fluorescence yield"
    ϕ_f::FT = 0
    "Non-photochemical quenching yeild"
    ϕ_n::FT = 0
    "Photochemical yield"
    ϕ_p::FT = 0
    "Dark adapted yield (`Kp=0`)"
    f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    f_m′::FT = 0
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    f_o′::FT = 0
    "Rate constant for thermal dissipation"
    k_d::FT = 0
    "Reversible NPQ rate constant (initially zero)"
    k_npq_rev::FT = 0
    "Rate constant for photochemistry"
    k_p::FT = 0
    "Non-Photochemical quenching "
    npq::FT = 0
    "Energy quenching"
    q_e::FT = 0
    "Photochemical quenching"
    q_p::FT = 0
    "max PSII yield (_k_npq_rev = 0, all RC open)"
    ϕ_psii_max::FT = 0

    # Cytochrome reaction center vars
    "Weight factor that PSI fluorescence reaches sensor (after reabsorption)"
    ϵ_1::FT = 0
    "Weight factor that PSII fluorescence reaches sensor (after reabsorption)"
    ϵ_2::FT = 1

    # cache variables
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-03: add C3VJP struct
#     2023-Oct-03: add C3Cyto struct
#     2023-Oct-03: add C4VJP struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the fields for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3VJP{FT}
    "State variables"
    state::C3VJPState{FT} = C3VJPState{FT}()
    "Auxilary variables"
    auxil::PSMAuxil{FT} = PSMAuxil{FT}()
end


"""

$(TYPEDEF)

Struct that contains the fields for C3 photosynthesis (Cytochrome model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3Cyto{FT}
    "State variables"
    state::C3CytoState{FT} = C3CytoState{FT}()
    "Auxilary variables"
    auxil::PSMAuxil{FT} = PSMAuxil{FT}()
end


"""

$(TYPEDEF)

Struct that contains the fields for C4 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4VJP{FT}
    "State variables"
    state::C4VJPState{FT} = C4VJPState{FT}()
    "Auxilary variables"
    auxil::PSMAuxil{FT} = PSMAuxil{FT}()
end

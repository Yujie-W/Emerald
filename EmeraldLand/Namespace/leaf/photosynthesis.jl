# This file contains the state and auxiliary variables for leaf photosynthesis

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-23: add C3CytoInfApTrait struct
#     2024-Jul-27: use modified TD for η_C and η_L
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis (Cytochrome model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3CytoInfApTrait{FT}
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
    TD_ηC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηCTDWang(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_L temperature dependency"
    TD_ηL::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηLTDWang(FT)

    # Constant coefficients
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

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add C3CytoTrait struct
#     2024-Jul-27: use modified TD for η_C and η_L
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis (Cytochrome model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3CytoTrait{FT}
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
    TD_ηC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηCTDWang(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_L temperature dependency"
    TD_ηL::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηLTDWang(FT)

    # Constant coefficients
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

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-22: add C3JBTrait struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis (Cytochrome model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3JBTrait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::MinimumColimit{FT} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::MinimumColimit{FT} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::SerialColimit{FT} = SerialColimit{FT}()

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

    # Constant coefficients
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

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add C3VJPTrait struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3VJPTrait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = ColimitJCLM(FT)

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
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.85
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for PSII photochemistry"
    K_PSII::FT = 4

    # Embedded structures
    "Fluorescence model"
    FLM::Union{KNFluoscenceModel{FT}, QLFluoscenceModel{FT}} = KNFluoscenceModel{FT}()

    # Prognostic variables
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
#     2024-Jul-22: add C3CLMTrait struct
#     2024-Jul-22: CLM colimitations are all quadratic
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3CLMTrait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::QuadraticColimit{FT} = ColimitCJCLMC3(FT)
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::QuadraticColimit{FT} = ColimitIPCLM(FT)
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::QuadraticColimit{FT} = ColimitJCLM(FT)

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
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.85
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for PSII photochemistry"
    K_PSII::FT = 4

    # Embedded structures
    "Fluorescence model"
    FLM::Union{KNFluoscenceModel{FT}, QLFluoscenceModel{FT}} = KNFluoscenceModel{FT}()

    # Prognostic variables
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
#     2024-Jul-22: add C3FvCBTrait struct
#     2024-Jul-22: FvCB colimitations are all minimum and quadratic (TODO: new algorithm required)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3FvCBTrait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::MinimumColimit{FT} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::MinimumColimit{FT} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::QuadraticColimit{FT} = QuadraticColimit{FT}(CURVATURE = 0.7)

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
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.85
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for PSII photochemistry"
    K_PSII::FT = 4

    # Embedded structures
    "Fluorescence model"
    FLM::Union{KNFluoscenceModel{FT}, QLFluoscenceModel{FT}} = KNFluoscenceModel{FT}()

    # Prognostic variables
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
#     2024-Apr-15: add C4CLMTrait struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C4 photosynthesis (CLM VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4CLMTrait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency"
    TD_KPEP::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = Q10TDKpepCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type  respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)

    # Constant coefficients
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.85
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for PSII photochemistry"
    K_PSII::FT = 4

    # Embedded structures
    "Fluorescence model"
    FLM::Union{KNFluoscenceModel{FT}, QLFluoscenceModel{FT}} = KNFluoscenceModel{FT}()

    # Prognostic variables
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add C4VJPTrait struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for C4 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4VJPTrait{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}, SquareColimit{FT}} = MinimumColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency"
    TD_KPEP::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KpepTDBoyd(FT)
    "[`AbstractTemperatureDependency`](@ref) type  respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vpmax temperature dependency"
    TD_VPMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VpmaxTDBoyd(FT)

    # Constant coefficients
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.85
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for PSII photochemistry"
    K_PSII::FT = 4

    # Embedded structures
    "Fluorescence model"
    FLM::Union{KNFluoscenceModel{FT}, QLFluoscenceModel{FT}} = KNFluoscenceModel{FT}()

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
    f_m::Vector{FT}
    "Light adapted yield (`Kp=0`)"
    f_m′::Vector{FT}
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::Vector{FT}
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

CanopyLayerPhotosystemAuxil(config::SPACConfiguration{FT}) where {FT} = CanopyLayerPhotosystemAuxil{FT}(
    a_c   = zeros(FT, config.DIM_PPAR_BINS+1),
    a_g   = zeros(FT, config.DIM_PPAR_BINS+1),
    a_i   = zeros(FT, config.DIM_PPAR_BINS+1),
    a_j   = zeros(FT, config.DIM_PPAR_BINS+1),
    a_n   = zeros(FT, config.DIM_PPAR_BINS+1),
    a_p   = zeros(FT, config.DIM_PPAR_BINS+1),
    e2c   = zeros(FT, config.DIM_PPAR_BINS+1),
    j     = zeros(FT, config.DIM_PPAR_BINS+1),
    j_pot = zeros(FT, config.DIM_PPAR_BINS+1),
    j_psi = zeros(FT, config.DIM_PPAR_BINS+1),
    η     = zeros(FT, config.DIM_PPAR_BINS+1),
    ϕ_f   = zeros(FT, config.DIM_PPAR_BINS+1),
    ϕ_p   = zeros(FT, config.DIM_PPAR_BINS+1),
    ϕ_f1  = zeros(FT, config.DIM_PPAR_BINS+1),
    ϕ_f2  = zeros(FT, config.DIM_PPAR_BINS+1),
    f_m   = zeros(FT, config.DIM_PPAR_BINS+1),
    f_m′  = zeros(FT, config.DIM_PPAR_BINS+1),
    f_o   = zeros(FT, config.DIM_PPAR_BINS+1),
    f_o′  = zeros(FT, config.DIM_PPAR_BINS+1),
    npq   = zeros(FT, config.DIM_PPAR_BINS+1),
    q_e   = zeros(FT, config.DIM_PPAR_BINS+1),
    q_p   = zeros(FT, config.DIM_PPAR_BINS+1),
    k_n   = zeros(FT, config.DIM_PPAR_BINS+1),
    k_p   = zeros(FT, config.DIM_PPAR_BINS+1)
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
    trait::Union{C3CytoInfApTrait{FT}, C3CytoTrait{FT}, C3CLMTrait{FT}, C3FvCBTrait{FT}, C3JBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}} = C3VJPTrait{FT}()
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
    trait::Union{C3CytoInfApTrait{FT}, C3CytoTrait{FT}, C3CLMTrait{FT}, C3FvCBTrait{FT}, C3JBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}} = C3VJPTrait{FT}()
    "State variables"
    state::Union{C3State{FT}, C4State{FT}} = C3State{FT}()
    "Auxilary variables"
    auxil::CanopyLayerPhotosystemAuxil{FT}
end;

CanopyLayerPhotosystem(config::SPACConfiguration{FT}) where {FT} = CanopyLayerPhotosystem{FT}(auxil = CanopyLayerPhotosystemAuxil(config));

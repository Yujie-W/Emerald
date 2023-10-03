# This file contains the state and auxilary variables for leaf photosynthesis

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
end

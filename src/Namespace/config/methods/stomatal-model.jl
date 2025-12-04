"""
Hierarchy of AbstractStomatalConductanceModel:
- AndereggSM
- BallBerrySM
- EllerSM
- GentineSM
- LeuningSM
- MedlynSM
- SperrySM
- WangSM
- Wang2SM

"""
abstract type AbstractStomatalConductanceModel{FT<:AbstractFloat} end;


"""
Struct for Anderegg stomatal model. The equation used for Anderegg type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{2aP + b}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
Base.@kwdef mutable struct AndereggSM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    A::FT = 0.5
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    B::FT = 2
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;


"""
Struct for Ball Berry stomatal model. The equation used for Ball-Berry type model is
```math
gs = g0 + g1 ⋅ RH ⋅ \\dfrac{A}{Cs}
```
"""
Base.@kwdef mutable struct BallBerrySM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT = 9
    "Beta function to force stomatal response to soil moisture"
    β::BetaFunction{FT} = BetaFunction{FT}()
    "Time constant for the prognostic stomatal conductance `[s]`"
    τ::FT = 600
end;


"""
Empty struct for Eller stomatal model. The equation used for Eller type model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
Base.@kwdef mutable struct EllerSM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;


"""
Struct for Gentine stomatal model. The equation used for Gentine type model is
```math
gs = g0 + g1 ⋅ \\dfrac{k_{leaf}}{k_{max}} ⋅ \\dfrac{A}{Ci}.
```
"""
Base.@kwdef mutable struct GentineSM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT = 9
    "Beta function to force stomatal response to soil moisture"
    β::BetaFunction{FT} = BetaFunction{FT}(FUNC = (x -> x), PARAM_X = BetaParameterKleaf(), PARAM_Y = BetaParameterG1())
    "Time constant for the prognostic stomatal conductance `[s]`"
    τ::FT = 600
end;


"""
Struct for Leuning stomatal model. The equation used for Leuning type model is
```math
gs = g0 + g1 ⋅ \\dfrac{A}{Cs - Γ^{*}} ⋅ \\dfrac{1}{1 + \\dfrac{VPD}{d0}}
```
"""
Base.@kwdef mutable struct LeuningSM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Fitting parameter of d/d0 below the fraction, same unit as vpd `[Pa]`"
    D0::FT = 3000
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT = 8
    "Beta function to force stomatal response to soil moisture"
    β::BetaFunction{FT} = BetaFunction{FT}()
    "Time constant for the prognostic stomatal conductance `[s]`"
    τ::FT = 600
end;


"""
Struct for Medlyn stomatal model. The equation used for Medlyn type model is
```math
gs = g0 + 1.6 ⋅ \\left( 1 + \\dfrac{g1}{\\sqrt{VPD}} \\right) ⋅ \\dfrac{A}{Ca}
```
"""
Base.@kwdef mutable struct MedlynSM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[sqrt(Pa)]`"
    G1::FT = 125
    "Beta function to force stomatal response to soil moisture"
    β::BetaFunction{FT} = BetaFunction{FT}()
    "Time constant for the prognostic stomatal conductance `[s]`"
    τ::FT = 600
end;


"""
Empty struct for Sperry stomatal model. The equation used for Sperry type model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A_{max}}{K_{max}}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
Base.@kwdef mutable struct SperrySM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;


"""
Empty struct for Wang stomatal model. The equation used for Wang type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{A}{E_{crit} - E}
```
"""
Base.@kwdef mutable struct WangSM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Fitness factor"
    F_FITNESS::FT = 0.2
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;


"""
Empty struct for a new Wang stomatal model modified from Anderegg model. The equation used for new Wang2SM type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{aAP}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
Base.@kwdef mutable struct Wang2SM{FT<:AbstractFloat} <: AbstractStomatalConductanceModel{FT}
    # General model information
    "Quadratic equation parameter `[MPa⁻²]`"
    A::FT = 0.1
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;

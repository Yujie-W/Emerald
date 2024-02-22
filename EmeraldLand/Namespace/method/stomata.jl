# This file contains the stomtal models related structs

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for stomatal models
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractStomataModel:
- [`AndereggSM`](@ref)
- [`BallBerrySM`](@ref)
- [`EllerSM`](@ref)
- [`GentineSM`](@ref)
- [`LeuningSM`](@ref)
- [`MedlynSM`](@ref)
- [`SperrySM`](@ref)
- [`WangSM`](@ref)
- [`Wang2SM`](@ref)

"""
abstract type AbstractStomataModel{FT<:AbstractFloat} end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Anderegg model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Anderegg stomatal model. The equation used for Anderegg type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{2aP + b}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AndereggSM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
    # General model information
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    A::FT = 0.5
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    B::FT = 2
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;

sync_state!(state_from::AndereggSM{FT}, state_to::AndereggSM{FT}) where {FT} = (
    state_to.A = state_from.A;
    state_to.B = state_from.B;
    state_to.K = state_from.K;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Ball Berry model
#     2022-Jul-07: add time constant
#     2022-Jul-20: rename Β and Τ (both greek) to β and τ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Ball Berry stomatal model. The equation used for Ball-Berry type model is
```math
gs = g0 + g1 ⋅ RH ⋅ \\dfrac{A}{Cs}
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BallBerrySM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
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

sync_state!(state_from::BallBerrySM{FT}, state_to::BallBerrySM{FT}) where {FT} = (
    state_to.G0 = state_from.G0;
    state_to.G1 = state_from.G1;
    sync_state!(state_from.β, state_to.β);
    state_to.τ = state_from.τ;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Eller model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Eller stomatal model. The equation used for Eller type model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
Base.@kwdef mutable struct EllerSM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
    # General model information
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;

sync_state!(state_from::EllerSM{FT}, state_to::EllerSM{FT}) where {FT} = (
    state_to.K = state_from.K;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Gentine model
#     2022-Jul-07: add field Β to generalize the model as others
#     2022-Jul-07: add time constant
#     2022-Jul-20: rename Β and Τ (both greek) to β and τ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Gentine stomatal model. The equation used for Gentine type model is
```math
gs = g0 + g1 ⋅ \\dfrac{k_{leaf}}{k_{max}} ⋅ \\dfrac{A}{Ci}.
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GentineSM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
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

sync_state!(state_from::GentineSM{FT}, state_to::GentineSM{FT}) where {FT} = (
    state_to.G0 = state_from.G0;
    state_to.G1 = state_from.G1;
    sync_state!(state_from.β, state_to.β);
    state_to.τ = state_from.τ;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Leuning model
#     2022-Jul-07: add time constant
#     2022-Jul-20: rename Β and Τ (both greek) to β and τ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Leuning stomatal model. The equation used for Leuning type model is
```math
gs = g0 + g1 ⋅ \\dfrac{A}{Cs - Γ^{*}} ⋅ \\dfrac{1}{1 + \\dfrac{VPD}{d0}}
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeuningSM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
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

sync_state!(state_from::LeuningSM{FT}, state_to::LeuningSM{FT}) where {FT} = (
    state_to.D0 = state_from.D0;
    state_to.G0 = state_from.G0;
    state_to.G1 = state_from.G1;
    sync_state!(state_from.β, state_to.β);
    state_to.τ = state_from.τ;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Medlyn model
#     2022-Jul-07: add time constant
#     2022-Jul-20: rename Β and Τ (both greek) to β and τ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Medlyn stomatal model. The equation used for Medlyn type model is
```math
gs = g0 + 1.6 ⋅ \\left( 1 + \\dfrac{g1}{\\sqrt{VPD}} \\right) ⋅ \\dfrac{A}{Ca}
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MedlynSM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
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

sync_state!(state_from::MedlynSM{FT}, state_to::MedlynSM{FT}) where {FT} = (
    state_to.G0 = state_from.G0;
    state_to.G1 = state_from.G1;
    sync_state!(state_from.β, state_to.β);
    state_to.τ = state_from.τ;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Sperry model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Sperry stomatal model. The equation used for Sperry type model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A_{max}}{K_{max}}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
Base.@kwdef mutable struct SperrySM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
    # General model information
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;

sync_state!(state_from::SperrySM{FT}, state_to::SperrySM{FT}) where {FT} = (
    state_to.K = state_from.K;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Wang model
#     2022-Jul-11: add fields: F_FITNESS, and ppar_mem
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Wang stomatal model. The equation used for Wang type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{A}{E_{crit} - E}
```
"""
Base.@kwdef mutable struct WangSM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
    # General model information
    "Fitness factor"
    F_FITNESS::FT = 0.2
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;

sync_state!(state_from::WangSM{FT}, state_to::WangSM{FT}) where {FT} = (
    state_to.F_FITNESS = state_from.F_FITNESS;
    state_to.K = state_from.K;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Wang model modified from Anderegg model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for a new Wang stomatal model modified from Anderegg model. The equation used for new Wang2SM type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{aAP}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Wang2SM{FT<:AbstractFloat} <: AbstractStomataModel{FT}
    # General model information
    "Quadratic equation parameter `[MPa⁻²]`"
    A::FT = 0.1
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end;

sync_state!(state_from::Wang2SM{FT}, state_to::Wang2SM{FT}) where {FT} = (
    state_to.A = state_from.A;
    state_to.K = state_from.K;

    return nothing
);

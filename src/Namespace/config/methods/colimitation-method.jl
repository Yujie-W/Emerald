"""
Hierarchy of `AbstractColimit`
- MinimumColimit
- QuadraticColimit
- SerialColimit

"""
abstract type AbstractColimit{FT<:AbstractFloat} end;


""" Empty structure to indicate minimum colimitation: `x = min(x₁, x₂)` """
struct MinimumColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end;


""" Structure to indicate quadratic colimitation (contains field `CURVATURE`): `θ⋅x² - (x₁ + x₂)⋅x + x₁x₂ = 0` """
Base.@kwdef mutable struct QuadraticColimit{FT<:AbstractFloat} <: AbstractColimit{FT}
    "Curvature factor"
    CURVATURE::FT = 0.98
end;


""" Empty structure to indicate serial colimitation: `x = 1 / (1/x₁ + 1/x₂)` """
struct SerialColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end;


""" Empty structure to indicate square colimitation: `x = x₁⋅x₂ / sqrt(x₁² + x₂²)` """
struct SquareColimit{FT<:AbstractFloat} <: AbstractColimit{FT} end;

# This file contains the parameters of the trace gasses

#######################################################################################################################################################################################################
#
# Changes to the types
# General:
#     2022-Oct-17: add trace gasses
#     2023-Jun-13: add trace gas struct for CH₄ and fields for N₂
#
#######################################################################################################################################################################################################
"""

`WaterPhysics` uses the multiple dispatch approach to calculate the temperature and pressure dependent physical properties of water and other molecules, such as CO₂. The trace molecules and mediums
    are catergorized to gas or liquid subject to a general type `AbstractTrace`. Hierarchy of `AbstractTrace`:
- [`AbstractTraceGas`](@ref)
- [`AbstractTraceLiquid`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTrace{FT<:AbstractFloat} end;


"""

The gas can be either the target trace molecule (e.g., when computing diffusive coefficient of CO₂ in water using [`diffusive_coefficient`](@ref)) or the medium (e.g., when computing diffusive
    coefficient or CO₂ in air using [`diffusive_coefficient`](@ref)). Currently, `WaterPhysics` supports the following subtypes of `AbstractTraceGas`:
- [`TraceGasAir`](@ref)
- [`TraceGasCH₄`](@ref)
- [`TraceGasCO₂`](@ref)
- [`TraceGasH₂O`](@ref)
- [`TraceGasO₂`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceGas{FT<:AbstractFloat} <: AbstractTrace{FT} end;


"""

The liquid can be either the medium for gas (e.g., when computing diffusive coefficient of CO₂ in water using [`diffusive_coefficient`](@ref)) or the target molecule (e.g., when computing surface
    tension of water using [`surface_tension`](@ref)). Currently. `WaterPhysics` supports the following subtypes of `AbstractTraceLiquid`:
- [`TraceLiquidH₂O`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceLiquid{FT<:AbstractFloat} <: AbstractTrace{FT} end;


"""

$(TYPEDEF)

Identity trace label for air.

"""
struct TraceGasAir{FT<:AbstractFloat} <: AbstractTraceGas{FT} end;


"""

$(TYPEDEF)

Identity trace label for gas phase CH₄.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceGasCH₄{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air::FT = 2.20e-5

    # related to diffusive coefficient in liquid water
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 1.49e-9
end;


"""

$(TYPEDEF)

Identity trace label for gas phase CO₂.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceGasCO₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air::FT = 2.82e-5 / 1.6

    # related to diffusive coefficient in liquid water
    "Hydrodynamic radius of the solute in `[m]`"
    a_298::FT = 1.68e-10
    "Coefficient to make temperature correction over ydrodynamic radius"
    a_a::FT = 2e-3
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 2.147813e-9
end;


"""

$(TYPEDEF)

Identity trace label for gas phase H₂O.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceGasH₂O{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air::FT = 2.82e-5
end;


"""

$(TYPEDEF)

Identity trace label for gas phase N₂.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceGasN₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air::FT = 2.12e-5

    # related to diffusive coefficient in liquid water
    "Hydrodynamic radius of the solute in `[m]`"
    a_298::FT = 1.90e-10
    "Coefficient to make temperature correction over ydrodynamic radius"
    a_a::FT = 2.2e-3
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 1.899062e-9
end;


"""

$(TYPEDEF)

Identity trace label for gas phase O₂.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceGasO₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air::FT = 1.76e-5

    # related to diffusive coefficient in liquid water
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 2.10e-9
end;


"""

$(TYPEDEF)

Identity trace label for liquid phase H₂O.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceLiquidH₂O{FT<:AbstractFloat} <: AbstractTraceLiquid{FT}
    # Surface tension related
    # γ = γ_k * (1-T/γ_T_c)^γ_exp * [1 - γ_cor*(1-T/γ_T_c)]
    "Surface tension coefficient correction"
    γ_cor::FT = 0.625
    "Surface tension coefficient exponent"
    γ_exp::FT = 1.256
    "Surface tension coefficient k in `[N m⁻¹]`"
    γ_k::FT = 0.2358
    "Surface tension at 298.15 K in `[N m⁻¹]`"
    γ_ref::FT = 0.07197220523
    "Surface tension coefficient T_crit in `[K]`"
    γ_T_c::FT = 647.096

    # Viscosity related
    # υ = A ⋅ exp( A/T + C⋅T + D⋅T^2 )
    "Viscosity coefficient A in `[Pa s]`"
    υ_A::FT = 1.856e-14
    "Viscosity coefficient B in `[K]`"
    υ_B::FT = 4209.0
    "Viscosity coefficient C in `[K⁻¹]`"
    υ_C::FT = 0.04527
    "Viscosity coefficient D in `[K⁻²]`"
    υ_D::FT = -3.376e-5
end;

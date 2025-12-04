"""
Hierarchy of AbstractTemperatureDependency:
- Arrhenius
- ArrheniusPeak
- ArrheniusPeak2
- Q10
- Q10Peak
- Q10PeakHT
- Q10PeakLTHT
"""
abstract type AbstractTemperatureDependency{FT<:AbstractFloat} end;


"""
An `Arrhenius` type struct using
```math
Y_1 = Y_0 \\cdot \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right)
```
"""
Base.@kwdef mutable struct Arrhenius{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
end;


"""
An `ArrheniusPeak` type struct using
```math
Y_1 = Y_0 \\cdot \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right)
          \\cdot \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```
"""
Base.@kwdef mutable struct ArrheniusPeak{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
    "Deactivation energy"
    ΔHD::FT
    "Entropy factor"
    ΔSV::FT
end;


"""
An `ArrheniusPeak2` type struct using
```math
Y_1 = Y_0 \\cdot \\min \\left(1, \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right) \\right)
          \\cdot \\min \\left(1, \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) } \\right)
```
"""
Base.@kwdef mutable struct ArrheniusPeak2{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
    "Deactivation energy"
    ΔHD::FT
    "Entropy factor"
    ΔSV::FT
end;


"""
A `Q10` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
```
"""
Base.@kwdef mutable struct Q10{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Power of Q10 correction"
    Q_10::FT
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
end;


"""
A `Q10Peak` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
          \\cdot \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```
"""
Base.@kwdef mutable struct Q10Peak{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Power of Q10 correction"
    Q_10::FT
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Deactivation energy"
    ΔHD::FT
    "Entropy factor"
    ΔSV::FT
end;


"""
A `Q10PeakHT` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
          \\cdot \\dfrac{ 1 }
                        { 1 + \\exp \\left( S_H * (T_1 - T_H) \\right) }
```
"""
Base.@kwdef mutable struct Q10PeakHT{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Power of Q10 correction"
    Q_10::FT
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Reference temperature to compute ΔT `[K]`"
    ΔT_REF::FT
    "Slope for ΔT `[K⁻¹]`"
    ΔT_SLOPE::FT
end;


"""
A `Q10PeakLTHT` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
          \\cdot \\dfrac{ 1 }{ 1 + \\exp \\left( S_H * (T_1 - T_H) \\right) }
          \\cdot \\dfrac{ 1 }{ 1 + \\exp \\left( S_L * (T_L - T_1) \\right) }
```
"""
Base.@kwdef mutable struct Q10PeakLTHT{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Power of Q10 correction"
    Q_10::FT
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Reference high temperature to compute ΔT `[K]`"
    ΔHT_REF::FT
    "Slope for high ΔT `[K⁻¹]`"
    ΔHT_SLOPE::FT
    "Reference low temperature to compute ΔT `[K]`"
    ΔLT_REF::FT
    "Slope for low ΔT `[K⁻¹]`"
    ΔLT_SLOPE::FT
end;

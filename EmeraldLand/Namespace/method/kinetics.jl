# This file contains the temperature dependency methods

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-13: migrate abstract temperature dependency type from Photosynthesis.jl
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractTemperatureDependency:
- [`Arrhenius`](@ref)
- [`ArrheniusPeak`](@ref)
- [`Q10`](@ref)

"""
abstract type AbstractTemperatureDependency{FT<:AbstractFloat} end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to Arrhenius
#     2022-Jan-13: define the struct mutable, use ΔHA directly in the struct, add field T_REF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An `Arrhenius` type struct using
```math
Y_1 = Y_0 \\cdot \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right)
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to ArrheniusPeak
#     2022-Jan-13: define the struct mutable, use ΔHA/ΔHD/ΔSV directly in the struct, add field T_REF/VAL_REF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An `ArrheniusPeak` type struct using
```math
Y_1 = Y_0 \\cdot \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right)
          \\cdot \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Oct-04: add struct ArrheniusPeak2
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An `ArrheniusPeak2` type struct using
```math
Y_1 = Y_0 \\cdot \\min \\left(1, \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right) \\right)
          \\cdot \\min \\left(1, \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) } \\right)
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to Q10
#     2022-Jan-14: make structure mutable
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

A `Q10` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jul-29: add struct Q10Peak
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

A `Q10Peak` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
          \\cdot \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Aug-01: add struct Q10PeakHT
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

A `Q10PeakHT` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
          \\cdot \\dfrac{ 1 }
                        { 1 + \\exp \\left( S_H * (T_1 - T_H) \\right) }
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Aug-01: add struct Q10PeakLTHT
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

A `Q10PeakLTHT` type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
          \\cdot \\dfrac{ 1 }{ 1 + \\exp \\left( S_H * (T_1 - T_H) \\right) }
          \\cdot \\dfrac{ 1 }{ 1 + \\exp \\left( S_L * (T_L - T_1) \\right) }
```

# Fields

$(TYPEDFIELDS)

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


#######################################################################################################################################################################################################
#
# Changes to the constructors
# General
#     2022-Jan-14: migrate from Photosynthesis.jl
#     2022-Feb-11: add temperature dependent Jmax and Vcmax TD from CLM
#     2024-Jul-27: add modified TD for η_C and η_L (value to be optimized)
# Sources
#     Lavigne and Ryan (1997) Growth and maintenance respiration rates of aspen, blackspruce and jack pine stems at northern and southern BOREAS sites
#     Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco‐limited photosynthesis
#     Boyd et al. (2001) Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
#     Leuning (2002) Temperature dependence of two parameters in a photosynthesis model
#     Kattge et al. (2007) Temperature acclimation in a biochemical model of photosynthesis: a reanalysis of data from 36 species
#     Sperry et al. (2019) The impact of rising CO2 and acclimation on the response of US forests to global warming
#     Johnson et al. (2021) The limiting factors and regulatory processes that control the environmental responses of C3, C3–C4 intermediate, and C4 photosynthesis
#     CLM5 Documentation. Chapter 9 Page 106
#
#######################################################################################################################################################################################################
KcTDBernacchi(FT)          = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 41.0264925, ΔHA = 79430.0);
KcTDCLM(FT)                = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 40.49     , ΔHA = 79430.0);
KoTDBernacchi(FT)          = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 28208.88  , ΔHA = 36380.0);
KoTDCLM(FT)                = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 27840.0   , ΔHA = 36380.0);
KpepTDBoyd(FT)             = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 16.0      , ΔHA = 36300.0);
KqTDJohnson(FT)            = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 300       , ΔHA = 37000.0);
RespirationTDBernacchi(FT) = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = NaN       , ΔHA = 46390.0);
VcmaxTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = NaN       , ΔHA = 65330.0);
VomaxTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = NaN       , ΔHA = 60110.0);
ΓStarTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 4.33164375, ΔHA = 37830.0);
ΓStarTDCLM(FT)             = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 4.275     , ΔHA = 37830.0);

JmaxTDBernacchi(FT)                 = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 57500.0, ΔHD = 439000.0, ΔSV = 1400.0);
JmaxTDCLM(FT, t::Number = T₂₅())    = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 50000.0, ΔHD = 200000.0, ΔSV = 659.70 - 0.75 * (t - T₀(FT)) );
JmaxTDLeuning(FT)                   = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 50300.0, ΔHD = 152044.0, ΔSV = 495.0 );
RespirationTDCLMC3(FT)              = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 46390.0, ΔHD = 150650.0, ΔSV = 490.0 );
VcmaxTDCLMC3(FT, t::Number = T₂₅()) = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 72000.0, ΔHD = 200000.0, ΔSV = 668.39 - 1.07 * (t - T₀(FT)) );
VcmaxTDLeuning(FT)                  = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 73637.0, ΔHD = 149252.0, ΔSV = 486.0 );
VpmaxTDBoyd(FT)                     = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 94800.0, ΔHD = 73300.0 , ΔSV = 250.0 );
ηCTDJohnson(FT)                     = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 1.0 , ΔHA = 0.0    , ΔHD = 220000.0, ΔSV = 710.0 );
ηLTDJohnson(FT)                     = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 0.75, ΔHA = 0.0    , ΔHD = 220000.0, ΔSV = 710.0 );
ηCTDWang(FT)                        = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 1.0 , ΔHA = 0.0    , ΔHD = 225100.0, ΔSV = 710.0 );
ηLTDWang(FT)                        = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 0.75, ΔHA = 0.0    , ΔHD = 225100.0, ΔSV = 710.0 );

Q10TDAngiosperm(FT) = Q10{FT}(Q_10 = 1.4, T_REF = T₂₅(FT), VAL_REF = 2 * 0.0140 / 8760 * 1000); # μmol CO2 mol⁻¹ C biomass s⁻¹
Q10TDGymnosperm(FT) = Q10{FT}(Q_10 = 1.7, T_REF = T₂₅(FT), VAL_REF = 2 * 0.0425 / 8760 * 1000); # μmol CO2 mol⁻¹ C biomass s⁻¹
Q10TDKpepCLM(FT)    = Q10{FT}(Q_10 = 2.0, T_REF = T₂₅(FT), VAL_REF = 0.2);

RespirationTDCLMC4(FT) = Q10PeakHT{FT}(Q_10 = 2.0, T_REF = T₂₅(FT), VAL_REF = NaN, ΔT_REF = 328.15, ΔT_SLOPE = 1.3);

VcmaxTDCLMC4(FT) = Q10PeakLTHT{FT}(Q_10 = 2.0, T_REF = T₂₅(FT), VAL_REF = NaN, ΔHT_REF = 313.15, ΔHT_SLOPE = 0.3, ΔLT_REF = 288.15, ΔLT_SLOPE = 0.2);


#######################################################################################################################################################################################################
#
# New parameters for the temperature dependency based on fitting A-Ci curves I collected
# TODO: make it default in the future after the paper is accepted
#
#######################################################################################################################################################################################################
ΓStarTDWang2024(FT) = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 4.56, ΔHA = 11800.0);

JmaxTDWang2024(FT, t::Number = T₂₅())  = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN   , ΔHA = 50000, ΔHD = 201000, ΔSV = 659.70 - 0.75 * (t - T₀(FT)));
KqTDWang2024(FT)                       = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 300   , ΔHA = 21900, ΔHD = 232000, ΔSV = 700);
VcmaxTDWang2024(FT, t::Number = T₂₅()) = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN   , ΔHA = 63000, ΔHD = 204000, ΔSV = 668.39 - 1.07 * (t - T₀(FT)));
ηCTDWang2024(FT)                       = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 2*3/14, ΔHA = 21900, ΔHD = 232000, ΔSV = 700);
ηLTDWang2024(FT)                       = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 3*3/14, ΔHA = 21900, ΔHD = 232000, ΔSV = 700);

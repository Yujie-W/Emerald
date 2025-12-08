"""
Abstract type for photosynthesis rate calculation methods for Ac:
- AcMethodC3VcmaxPi
- AcMethodC4Vcmax
"""
abstract type AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac):
```math
A_\\text{c} = V_{\\text{cmax}} \\dfrac{(P_{\\text{i}} - \\Gamma)}{(P_{\\text{i}} + K_{\\text{m}})}
```
"""
struct AcMethodC3VcmaxPi <: AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac):
```math
A_\\text{c} = V_{\\text{cmax}}
```
"""
struct AcMethodC4Vcmax <: AbstractAcMethod end;


"""
Abstract type for photosynthesis rate calculation methods for Aj:
- AjMethodC3JmaxPi
- AjMethodC3VqmaxPi
- AjMethodC4JPSII
"""
abstract type AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj):
```math
A_\\text{j} = \\dfrac{J}{4} \\dfrac{(P_{\\text{i}} - \\Gamma)}{(P_{\\text{i}} + 2 \\Gamma)}
```
"""
struct AjMethodC3JmaxPi <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:
```math
J = \\dfrac{J_{\\text{PSI}}}{\\eta}
A_{\\text{j}} = \\dfrac{J}{4} \\dfrac{(P_{\\text{i}} - \\Gamma)}{(P_{\\text{i}} + 2 \\Gamma)}
```
"""
struct AjMethodC3VqmaxPi <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:
```math
Aj = \\dfrac{J}{6}
```
"""
struct AjMethodC4JPSII <: AbstractAjMethod end;


"""
Abstract type for photosynthesis rate calculation methods for Ap:
- ApMethodC3Inf
- ApMethodC3Vcmax
- ApMethodC4VcmaxPi
- ApMethodC4VpmaxPi
"""
abstract type AbstractApMethod end;


"""
Method for models without TPU limitation:
```math
A_{\\text{p}} = \\infty
```
"""
struct ApMethodC3Inf <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:
```math
A_{\\text{p}} = \\dfrac{V_{\\text{cmax}}}{2}
```
"""
struct ApMethodC3Vcmax <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:
```math
A_{\\text{p}} = k \\cdot V_{\\text{cmax}} \\cdot P_{\\text{i}}
```
"""
struct ApMethodC4VcmaxPi <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vpmax:
```math
A_{\\text{p}} = V_{\\text{pmax}} \\dfrac{P_{\\text{i}}}{(P_{\\text{i}} + K_{\\text{pep}})}
```
"""
struct ApMethodC4VpmaxPi <: AbstractApMethod end;

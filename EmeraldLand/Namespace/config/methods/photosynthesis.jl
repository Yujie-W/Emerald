"""
Abstract type for photosynthesis rate calculation methods for Ac:
- AcMethodC3VcmaxPi
- AcMethodC4Vcmax
"""
abstract type AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac):

    Ac = Vcmax * (Pi - Γ) / (Pi + Km)

"""
struct AcMethodC3VcmaxPi <: AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac):

    Ac = Vcmax

"""
struct AcMethodC4Vcmax <: AbstractAcMethod end;


"""
Abstract type for photosynthesis rate calculation methods for Aj:
- AjMethodC3JmaxPi
- AjMethodC3VqmaxPi
- AjMethodC4JPSII
"""
abstract type AbstractAjMethod <: AbstractAcMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj):

    Aj = J/4 * (Pi - Γ) / (Pi + 2 * Γ)

"""
struct AjMethodC3JmaxPi <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:

    J = J_PSI / η
    Aj = J / 4 * (Pi - Γ) / (Pi + 2 * Γ)

"""
struct AjMethodC3VqmaxPi <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:

    Aj = J / 6

"""
struct AjMethodC4JPSII <: AbstractAjMethod end;


"""
Abstract type for photosynthesis rate calculation methods for Ap:
- ApMethodC3Inf
- ApMethodC3Vcmax
- ApMethodC4VcmaxPi
- ApMethodC4VpmaxPi
"""
abstract type AbstractApMethod <: AbstractAjMethod end;


"""
Method for models without TPU limitation:

    Ap = Inf

"""
struct ApMethodC3Inf <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:

    Ap = Vcmax / 2

"""
struct ApMethodC3Vcmax <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:

    Ap = k * Vcmax * Pi

"""
struct ApMethodC4VcmaxPi <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vpmax:

    Ap = Vpmax * Pi / (Pi + Kpep)

"""
struct ApMethodC4VpmaxPi <: AbstractApMethod end;

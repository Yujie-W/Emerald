# This file contains the methods for computing Ac, Aj, and Ap

#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2024-Jul-31: add method structs for Ac, Aj, and Ap
#                  AcMethodC3VcmaxPi, AcMethodC4Vcmax
#                  AjMethodC3JmaxPi, AjMethodC3VqmaxPi, AjMethodC4JPSII
#                  ApMethodC3Inf, ApMethodC3Vcmax, ApMethodC4Vcmax, ApMethodC4Vpmax
#
#######################################################################################################################################################################################################
"""
Method to compute rubisco-limited photosynthesis rate (Ac):

    Ac = Vcmax * (Pi - Γ) / (Pi + Km)

"""
struct AcMethodC3VcmaxPi end;

"""
Method to compute rubisco-limited photosynthesis rate (Ac):

    Ac = Vcmax

"""
struct AcMethodC4Vcmax end;

"""
Method to compute electron transport-limited photosynthesis rate (Aj):

    Aj = J/4 * (Pi - Γ) / (Pi + 2 * Γ)

"""
struct AjMethodC3JmaxPi end;

"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:

    J = J_PSI / η
    Aj = J / 4 * (Pi - Γ) / (Pi + 2 * Γ)

"""
struct AjMethodC3VqmaxPi end;

"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:

    Aj = J / 6

"""
struct AjMethodC4JPSII end;

"""
Method for models without TPU limitation:

    Ap = Inf

"""
struct ApMethodC3Inf end;

"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:

    Ap = Vcmax / 2

"""
struct ApMethodC3Vcmax end;

"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:

    Ap = k * Vcmax * Pi

    """
struct ApMethodC4Vcmax end;

"""
Method to compute product-limited photosynthesis rate (Ap) from Vpmax:

    Ap = Vpmax * Pi / (Pi + Kpep)

"""
struct ApMethodC4Vpmax end;

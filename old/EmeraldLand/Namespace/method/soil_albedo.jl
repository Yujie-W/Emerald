# This file contains methods to compute soil albedo

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2023-Oct-26: add abstract type for soil albedo
#     2024-Apr-19: add new method to prescribe soil albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilAlbedo:
- [`SoilAlbedoPrescribe`](@ref)
- [`SoilAlebedoBroadbandCLM`](@ref)
- [`SoilAlebedoBroadbandCLIMA`](@ref)
- [`SoilALbedoHyperspectralCLM`](@ref)
- [`SoilALbedoHyperspectralCLIMA`](@ref)

"""
abstract type AbstractSoilAlbedo end;


""" Broadband soil albedo algorithm from CLM """
struct SoilAlbedoPrescribe <: AbstractSoilAlbedo end;


""" Broadband soil albedo algorithm from CLM """
struct SoilAlbedoBroadbandCLM <: AbstractSoilAlbedo end;


""" Broadband soil albedo algorithm from CLIMA """
struct SoilAlbedoBroadbandCLIMA <: AbstractSoilAlbedo end;


""" Hyperspectral soil albedo fitted from CLM broadband soil albedo """
struct SoilAlbedoHyperspectralCLM <: AbstractSoilAlbedo end;


""" Hyperspectral soil albedo fitted from CLIMA broadband soil albedo """
struct SoilAlbedoHyperspectralCLIMA <: AbstractSoilAlbedo end;

"""
Hierarchy of AbstractSoilAlbedo:
- SoilAlbedoPrescribe
- SoilAlebedoBroadbandCLM
- SoilAlebedoBroadbandCLIMA
- SoilALbedoHyperspectralCLM
- SoilALbedoHyperspectralCLIMA
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

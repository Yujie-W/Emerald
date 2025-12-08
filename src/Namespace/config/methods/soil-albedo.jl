"""
Hierarchy of AbstractSoilAlbedo:
- SoilAlbedoPrescribe
- SoilAlebedoBroadbandCLM
- SoilAlebedoBroadbandCLIMA
- SoilALbedoHyperspectralCLM
- SoilALbedoHyperspectralCLIMA
- SoilAlbedoHyperspectralAsh
"""
abstract type AbstractSoilAlbedo end;


""" Broadband soil albedo algorithm from CLM """
struct SoilAlbedoPrescribe <: AbstractSoilAlbedo end;


""" Broadband soil albedo algorithm from CLM """
struct SoilAlbedoBroadbandCLM <: AbstractSoilAlbedo end;


""" Broadband soil albedo algorithm from CLIMA """
struct SoilAlbedoBroadbandCLIMA <: AbstractSoilAlbedo end;

""" Hyperspectral soil albedo method based on soil water content and ash coverage """
struct SoilAlbedoHyperspectralAsh <: AbstractSoilAlbedo end;


""" Hyperspectral soil albedo fitted from CLM broadband soil albedo """
struct SoilAlbedoHyperspectralCLM <: AbstractSoilAlbedo end;


""" Hyperspectral soil albedo fitted from CLIMA broadband soil albedo """
struct SoilAlbedoHyperspectralCLIMA <: AbstractSoilAlbedo end;

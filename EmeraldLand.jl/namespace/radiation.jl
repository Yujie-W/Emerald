#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-10: refactor the structure with renamed fields
#     2021-Aug-10: add a constructor within the structure to avoid external initialization
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#     2023-Jun-16: remove fields of DIM_*
#     2023-Sep-11: add field ΔΛ_SIF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct WaveLengthSet{FT<:AbstractFloat}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # Constants
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT} = FT[700, 2500]
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT} = FT[400, 750]
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT} = FT[640, 850]
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT} = FT[400, 750]
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT} = read_nc(DATASET, "WL")
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::Vector{FT} = read_nc(DATASET, "WL_LOWER")
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::Vector{FT} = read_nc(DATASET, "WL_UPPER")

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int} = findall( WL_NIR[1] .<= Λ .<= WL_NIR[2] )
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int} = findall( WL_PAR[1] .<= Λ .<= WL_PAR[2] )
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int} = findall( WL_SIF[1] .<= Λ .<= WL_SIF[2] )
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int} = findall( WL_SIFE[1] .<= Λ .<= WL_SIFE[2] )

    # Constants based on the ones above
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT} = Λ_UPPER .- Λ_LOWER
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for SIF `[nm]`"
    ΔΛ_SIF::Vector{FT} = ΔΛ[IΛ_SIF]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT} = Λ[IΛ_PAR]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT} = Λ[IΛ_SIFE]
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for incoming solar radiation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractRadiation:
- [`BroadbandRadiation`](@ref)
- [`HyperspectralRadiation`](@ref)

"""
abstract type AbstractRadiation{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-15: add broadband solar radiation
#     2022-Jul-19: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores broadband radiation information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandRadiation{FT<:AbstractFloat} <: AbstractRadiation{FT}
    # prognostic variables that change with time
    "Diffuse radiation from NIR region `[W m⁻²]`"
    e_diffuse_nir::FT = 0
    "Diffuse radiation from PAR region `[W m⁻²]`"
    e_diffuse_par::FT = 0
    "Direct radiation from NIR region `[W m⁻²]`"
    e_direct_nir::FT = 0
    "Direct radiation from PAR region `[W m⁻²]`"
    e_direct_par::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Oct-22: refactor the structure with renamed fields
#     2021-Oct-22: add a constructor to define the structure from wavelength sets and prescribed wave shape
#     2022-Jun-15: change the order of variables
#     2022-Jul-20: use kwdef for the constructor
#     2023-Sep-14: add constructor to define the structure from a NetCDF dataset
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralRadiation{FT<:AbstractFloat} <: AbstractRadiation{FT}
    # Prognostic variables
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    e_diffuse::Vector{FT}
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    e_direct::Vector{FT}
end

HyperspectralRadiation{FT}(dataset::String = LAND_2021) where {FT} = HyperspectralRadiation{FT}(read_nc(dataset, "E_DIFF"), read_nc(dataset, "E_DIR"));

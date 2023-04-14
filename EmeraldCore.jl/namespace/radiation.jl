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
abstract type AbstractRadiation{FT} end


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
Base.@kwdef mutable struct BroadbandRadiation{FT} <: AbstractRadiation{FT}
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
#     2022-Apr-14: add DIM_WL to radiation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralRadiation{FT,DIM_WL} <: AbstractRadiation{FT}
    # Prognostic variables
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    e_diffuse::Vector{FT} = zeros(FT,DIM_WL)
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    e_direct::Vector{FT} = zeros(FT,DIM_WL)
end

HyperspectralRadiation{FT,DIM_WL}(dset::String) where {FT,DIM_WL} = HyperspectralRadiation{FT,DIM_WL}(read_nc(dset, "E_DIFF"), read_nc(dset, "E_DIR"));

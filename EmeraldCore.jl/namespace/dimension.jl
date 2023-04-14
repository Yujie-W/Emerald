#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-14: move all dimensions to SPACDimension
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Global configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct SPACDimension
    "Dimension of air layers"
    DIM_AIR::Int = 25
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of canopy layers"
    DIM_CANOPY::Int = 20
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Number of wavelength bins for NIR"
    DIM_NIR::Int
    "Number of wavelength bins for PAR"
    DIM_PAR::Int
    "Dimension of root layers"
    DIM_ROOT::Int = 4
    "Dimension of SIF wave length bins"
    DIM_SIF::Int
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int
    "Dimension of soil layers"
    DIM_SOIL::Int = 4
    "Dimension of short wave length bins"
    DIM_WL::Int
    "Dimension of xylem elements"
    DIM_XYLEM::Int = 5
end

SPACDimension(dset::String, wl_nir::Vector{FT}, wl_par::Vector{FT}, wl_sif::Vector{FT}, wl_sife::Vector{FT}) where {FT} = (
    _λ = read_nc(dset, "WL");

    return SPACDimension(
                DIM_NIR  = length(findall(wl_nir[1]  .<= _λ .<= wl_nir[2])),
                DIM_PAR  = length(findall(wl_par[1]  .<= _λ .<= wl_par[2])),
                DIM_SIF  = length(findall(wl_sif[1]  .<= _λ .<= wl_sif[2])),
                DIM_SIFE = length(findall(wl_sife[1] .<= _λ .<= wl_sife[2])),
                DIM_WL   = length(_λ),
    )
);

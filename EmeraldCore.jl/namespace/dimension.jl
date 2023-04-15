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
    "Dimension of temperature history"
    DIM_HISTORY::Int = 240
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

SPACDimension(gcf::GeneralConfiguration, ncanopy::Int = 20) = (
    _λ = read_nc(gcf.DATASET, "WL");

    return SPACDimension(
                DIM_CANOPY = ncanopy,
                DIM_NIR    = length(findall(gcf.WL_NIR[1]  .<= _λ .<= gcf.WL_NIR[2])),
                DIM_PAR    = length(findall(gcf.WL_PAR[1]  .<= _λ .<= gcf.WL_PAR[2])),
                DIM_SIF    = length(findall(gcf.WL_SIF[1]  .<= _λ .<= gcf.WL_SIF[2])),
                DIM_SIFE   = length(findall(gcf.WL_SIFE[1] .<= _λ .<= gcf.WL_SIFE[2])),
                DIM_WL     = length(_λ),
    );
);

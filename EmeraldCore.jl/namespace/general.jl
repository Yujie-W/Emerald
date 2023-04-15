#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Apr-14: make a general struct that will not be used in any SPAC calculations
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

General configuration of SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GeneralConfiguration
    # Constants
    "Wavelength limits for NIR used for soil albedo `[nm]`"
    WL_NIR::Tuple{Number,Number} = (700, 2500)
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Tuple{Number,Number} = (400, 750)
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Tuple{Number,Number} = (640, 850)
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Tuple{Number,Number} = (400, 750)

    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021
end

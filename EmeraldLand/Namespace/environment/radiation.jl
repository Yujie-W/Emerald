# This file contains the struct for solar radiation profile

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Oct-22: refactor the structure with renamed fields
#     2021-Oct-22: add a constructor to define the structure from wavelength sets and prescribed wave shape
#     2023-Sep-14: add constructor to define the structure from a NetCDF dataset
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields

$(TYPEDFIELDS)

"""
mutable struct ShortwaveRadiation{FT<:AbstractFloat}
    # Prognostic variables
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    e_diffuse::Vector{FT}
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    e_direct::Vector{FT}
end;

ShortwaveRadiation(config::SPACConfiguration{FT}) where {FT} = ShortwaveRadiation{FT}(read_nc(config.DATASET, "E_DIFF"), read_nc(config.DATASET, "E_DIR"));

broadband_radiation(FT) = ShortwaveRadiation{FT}(zeros(FT,2), zeros(FT,2));

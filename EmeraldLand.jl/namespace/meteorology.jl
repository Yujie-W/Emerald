#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jul-14: add Meteorology struct to store meteorological data
#     2022-Jul-14: add field t_precip
#     2022-Jul-20: add field t_air
#     2023-Apr-13: move RAD_LW and RAD_SW from MultiLayerSPAC
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores meteorological information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Meteorology{FT<:AbstractFloat}
    # Prognostic variables
    "Downwelling longwave radiation `[W m⁻²]`"
    rad_lw::FT = 100
    "Downwelling shortwave radiation"
    rad_sw::HyperspectralRadiation{FT} = HyperspectralRadiation{FT}()
    "Precipitation in form of rain (before interception) `[mol m⁻²]`"
    rain::FT = 0
    "Precipitation in form of snow (before interception) `[mol m⁻²]`"
    snow::FT = 0
    "Air temperature as the boundary condition for canopy airspace `[K]`"
    t_air::FT = T₂₅()
    "Precipitation temperature `[K]`"
    t_precip::FT = T₂₅()
end

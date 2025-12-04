#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-31: add PhotosystemsRateConstants
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores rate constants for photosystems (do not split PSI and PSII)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct PhotosystemsRateConstants{FT}
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.85
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for PSII photochemistry"
    K_P::FT = 4
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-31: add PhotosystemIRateConstants
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores rate constants for photosystem I

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct PhotosystemIRateConstants{FT}
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.55
    "Rate constant of fluorescence `[ns⁻¹]`"
    K_F::FT = 0.05
    "Rate constant of photochemistry for PS I `[ns⁻¹]`"
    K_P::FT = 14.5
    "Rate constant of regulated heat loss via oxidized PS I center `[s⁻¹]`"
    K_X::FT = 14.5
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Jul-31: add PhotosystemIIRateConstants
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores rate constants for photosystem II

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct PhotosystemIIRateConstants{FT}
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.55
    "Rate constant of fluorescence `[ns⁻¹]`"
    K_F::FT = 0.05
    "Rate constant of photochemistry for PS II `[ns⁻¹]`"
    K_P::FT = 4.5
    "Rate constant of excitation sharing for PS II `[ns⁻¹]`"
    K_U::FT = 2
end;

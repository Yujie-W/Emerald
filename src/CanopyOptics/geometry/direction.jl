# this file contains functions to compute the fraction of backward and forward light directions

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Mar-01: add function f_adaxial to compute the fraction of light from adaxial side (for both direct and diffuse)
#
#######################################################################################################################################################################################################
"""

    f_adaxial(θ_sza::FT, θ_incl::FT) where {FT}
    f_adaxial(θ_incl::FT) where {FT}

Compute the fraction of light from adaxial side (if θ_sza is not provided, compute that for isotropic light), given
- `θ_sza` solar zenith angle
- `θ_incl` leaf inclination angle

"""
function f_adaxial end;

f_adaxial(θ_sza::FT, θ_incl::FT) where {FT} = (
    # if sza = 90, half light is from adaxial
    if θ_sza == 90
        return FT(0.5)
    end;

    # if sza = 0 and incl = 90, half light is from adaxial
    if θ_sza == 0 || θ_incl == 90
        return FT(0.5)
    end;

    # if sza <= 90 - incl, all light is from adaxial
    if θ_sza + θ_incl <= 90
        return FT(1)
    end;

    # otherwise, the fraction is computed based on the geometry
    return 1 - acosd(1 / (tand(θ_sza) * tand(θ_incl))) / 180
);

f_adaxial(θ_incl::FT) where {FT} = (
    # compute the mean extinction coefficient for diffuse solar radiation from 18 angles
    Σf::FT = 0;
    Σs::FT = 0;
    for θ_sza in 0:1:89.9
        Σf += f_adaxial(θ_sza + FT(0.5), θ_incl) * sind(θ_sza + FT(0.5));
        Σs += sind(θ_sza + FT(0.5));
    end;

    return Σf / Σs
);

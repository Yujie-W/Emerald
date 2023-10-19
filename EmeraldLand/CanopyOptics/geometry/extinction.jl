# This file contains function to compute the extinction coefficients

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-07: add relative azimuth angle control
#     2022-Jun-16: add method for direct radiation
#     2022-Jun-16: add method for diffuse radiation
#     2022-Jun-21: add a sin(sza) weight to the diffuse light coefficient
# Bug fixes
#     2023-Mar-16: adjust _sb and _sf if negative value appears
# Sources
#     Verhoef (1998) Theory of radiative transfer models applied in optical remote sensing of vegetation canopies. Chapter 7
#
#######################################################################################################################################################################################################
"""

    extinction_coefficient(sza::FT, lia::FT) where {FT}
    extinction_coefficient(lia::FT) where {FT}

Return the extinction coefficient for direct radiation (when sza provided) and diffuse radiation (sza not provided), given
- `sza` Solar zenith angle in `°`
- `lia` Leaf inclination angle in `°`

"""
function extinction_coefficient end;

extinction_coefficient(sza::FT, lia::FT) where {FT} = (
    Cs = cosd(lia) * cosd(sza);
    Ss = sind(lia) * sind(sza);
    βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));

    return 2 / FT(π) / cosd(sza) * (Cs * (βs - FT(π) / 2) + Ss * sin(βs))
);

extinction_coefficient(lia::FT) where {FT} = (
    # compute the mean extinction coefficient for diffuse solar radiation from 18 angles
    kd::FT = 0;
    for sza in 0:5:89
        kd += extinction_coefficient(sza + FT(2.5), lia) * sind(sza);
    end;

    return kd / 18
);

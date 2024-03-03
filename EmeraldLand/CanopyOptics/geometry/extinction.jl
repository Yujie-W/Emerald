# This file contains function to compute the extinction coefficients

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-07: add relative azimuth angle control
#     2022-Jun-16: add method for direct radiation
#     2022-Jun-16: add method for diffuse radiation
#     2022-Jun-21: add a sin(θ_sza) weight to the diffuse light coefficient
# Bug fixes
#     2023-Mar-16: adjust _sb and _sf if negative value appears
#     2024-Mar-01: weight the coefficient based on sind(θ_sza+2.5) rather than 1 (non-breaking because this method has not been used anywhere)
# Sources
#     Verhoef (1998) Theory of radiative transfer models applied in optical remote sensing of vegetation canopies. Chapter 7
#
#######################################################################################################################################################################################################
"""

    extinction_coefficient(θ_sza::FT, θ_incl::FT) where {FT}
    extinction_coefficient(θ_incl::FT) where {FT}

Return the extinction coefficient for direct radiation (when θ_sza provided) and diffuse radiation (θ_sza not provided), given
- `θ_sza` Solar zenith angle in `°`
- `θ_incl` Leaf inclination angle in `°`

"""
function extinction_coefficient end;

extinction_coefficient(θ_sza::FT, θ_incl::FT) where {FT} = (
    Cs = cosd(θ_incl) * cosd(θ_sza);
    Ss = sind(θ_incl) * sind(θ_sza);
    βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));

    return 2 / FT(π) / cosd(θ_sza) * (Cs * (βs - FT(π) / 2) + Ss * sin(βs))
);

extinction_coefficient(θ_incl::FT) where {FT} = (
    # compute the mean extinction coefficient for diffuse solar radiation from 18 angles
    Σk::FT = 0;
    Σs::FT = 0;
    for θ_sza in 0:1:89.9
        Σk += extinction_coefficient(θ_sza + FT(0.5), θ_incl) * sind(θ_sza + FT(0.5));
        Σs += sind(θ_sza + FT(0.5));
    end;

    return Σk / Σs
);

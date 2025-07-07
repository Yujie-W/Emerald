# This file contains function to compute the extinction coefficients

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-07: add relative azimuth angle control
#     2022-Jun-16: add method for direct radiation
#     2022-Jun-16: add method for diffuse radiation
#     2022-Jun-21: add a sin(θ_za) weight to the diffuse light coefficient
#     2024-Sep-04: add method to compute extinction coefficient with clumping index effect
#     2024-Oct-16: do not weigh the extinction coefficient for diffuse radiation here, loop through it in the canopy_structure_aux! and canopy_structure! functions
# Bug fixes
#     2023-Mar-16: adjust _sb and _sf if negative value appears
#     2024-Mar-01: weight the coefficient based on sind(θ_za+2.5) rather than 1 (non-breaking because this method has not been used anywhere)
# Sources
#     Verhoef (1998) Theory of radiative transfer models applied in optical remote sensing of vegetation canopies. Chapter 7
#
#######################################################################################################################################################################################################
"""

    extinction_coefficient(θ_za::FT, θ_incl::FT) where {FT}
    extinction_coefficient(θ_incl::FT) where {FT}

Return the extinction coefficient for direct radiation, given
- `θ_za` (Solar ot Viewer) zenith angle in `°`
- `θ_incl` Leaf inclination angle in `°`

"""
function extinction_coefficient end;

extinction_coefficient(θ_za::FT, θ_incl::FT) where {FT} = (
    Cs = cosd(θ_incl) * cosd(θ_za);
    Ss = sind(θ_incl) * sind(θ_za);
    βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));

    return 2 / FT(π) / cosd(θ_za) * (Cs * (βs - FT(π) / 2) + Ss * sin(βs))
);

extinction_coefficient(θ_za::FT, θ_incl::FT, ci::FT) where {FT} = ci * extinction_coefficient(θ_za, θ_incl);

extinction_coefficient(θ_za::FT, θ_incl::FT, ci::ClumpingIndex{FT}) where {FT} =  extinction_coefficient(θ_za, θ_incl, ci.ci_0 * (1 - ci.ci_1 * cosd(θ_za + FT(0.5))));

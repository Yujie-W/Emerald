# This file contains functions to compute ∂T∂E for nocturnal stomatal opening

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂T∂E
#     2023-Oct-16: use effective emissivity in the calculation
#
#######################################################################################################################################################################################################
"""

    ∂T∂E(eff_ϵ::FT, t::FT, width::FT, wind::FT) where {FT}

Returns the marginal increase in leaf temperature per transpiration rate (per canopy layer), given
- `eff_ϵ` Effective emissivity used to compute the longwave radiation emission (effective LAI * leaf emissivity)
- `t` Leaf temperature
- `width` Leaf width
- `wind` Wind speed

"""
function ∂T∂E(eff_ϵ::FT, t::FT, width::FT, wind::FT) where {FT}
    λ = latent_heat_vapor(t) * M_H₂O(FT);
    g = FT(1.4) * FT(0.135) * sqrt(wind / (FT(0.72) * width));
    d = 2 * CP_D_MOL(FT) * g + 2 * eff_ϵ * 4 * K_STEFAN(FT) * t ^ 3;

    return λ / d
end;

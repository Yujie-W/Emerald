# This file contains functions to compute ∂R∂E for nocturnal stomatal opening

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂R∂E
#     2023-Oct-16: use effective emissivity in the calculation
#
#######################################################################################################################################################################################################
"""

    ∂R∂E(leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT}

Returns the marginal increase in leaf respiration rate per transpiration rate (per leaf area), given
- `leaf` `Leaf`, `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `eff_ϵ` Effective emissivity used to compute the longwave radiation emission (effective LAI * leaf emissivity)

"""
function ∂R∂E(leaf::Leaf{FT}, air::AirLayer{FT}, eff_ϵ::FT) where {FT}
    return ∂R∂T(leaf) * ∂T∂E(eff_ϵ, leaf.energy.auxil.t, leaf.bio.state.width, air.auxil.wind) * leaf.xylem.state.area
end;

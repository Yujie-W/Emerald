# This file contains function related to the soil water infiltration

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-06: add function soil water infiltration
#
#######################################################################################################################################################################################################
"""

    soil_water_infiltration!(spac::MultiLayerSPAC{FT}) where {FT}

Update the marginal increase of soil water content per layer, given
- `spac` `MultiLayerSPAC` SPAC

"""
function soil_water_infiltration!(spac::MultiLayerSPAC{FT}) where {FT}
    (; METEO, SOIL_BULK, SOILS) = spac;

    # update k, δψ, and flow rate among layers
    SOILS[1].auxil.∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[1].auxil.δz;

    N = length(SOIL_BULK.auxil.q);
    for i in 1:N
        SOILS[i  ].auxil.∂θ∂t -= SOIL_BULK.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[i  ].auxil.δz;
        SOILS[i+1].auxil.∂θ∂t += SOIL_BULK.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[i+1].auxil.δz;
    end;

    return nothing
end;

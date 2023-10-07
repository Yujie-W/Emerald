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

    # update the soil bulk auxiliary variables related to soil water and enrgy transfer between adjacent layers
    N = length(SOIL_BULK.auxil.k);
    for i in 1:N
        SOIL_BULK.auxil.k[i]        = 1 / (2 / SOILS[i].auxil.k + 2 / SOILS[i+1].auxil.k);
        SOIL_BULK.auxil.δψ[i]       = SOILS[i].auxil.ψ - SOILS[i+1].auxil.ψ + ρg_MPa(FT) * (SOILS[i].auxil.z - SOILS[i+1].auxil.z);
        SOIL_BULK.auxil.q[i]        = SOIL_BULK.auxil.k[i] * SOIL_BULK.auxil.δψ[i];
        SOIL_BULK.auxil.λ_layers[i] = 1 / (2 / SOILS[i].auxil.λ_soil_water + 2 / SOILS[i+1].auxil.λ_soil_water);
        SOIL_BULK.auxil.δt[i]       = SOILS[i].auxil.t - SOILS[i+1].auxil.t;
        SOIL_BULK.auxil.q_layers[i] = SOIL_BULK.auxil.λ_layers[i] * SOIL_BULK.auxil.δt[i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (SOIL_BULK.auxil.q[i] > 0) && (SOILS[i+1].state.θ >= SOILS[i+1].state.vc.Θ_SAT)
            SOIL_BULK.auxil.q[i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (SOIL_BULK.auxil.q[i] < 0) && (SOILS[i].state.θ >= SOILS[i].state.vc.Θ_SAT)
            SOIL_BULK.auxil.q[i] = 0;
        end;

        # if both layers are oversaturated, move the oversaturated part from lower layer to upper layer
        if (SOILS[i].state.θ >= SOILS[i].state.vc.Θ_SAT) && (SOILS[i+1].state.θ > SOILS[i+1].state.vc.Θ_SAT)
            SOIL_BULK.auxil.q[i] = -1 * (SOILS[i+1].state.θ - SOILS[i+1].state.vc.Θ_SAT) * SOILS[i+1].auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);
        end;
    end;

    # update k, δψ, and flow rate among layers
    SOILS[1].auxil.∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[1].auxil.δz;

    N = length(SOIL_BULK.auxil.q);
    for i in 1:N
        SOILS[i  ].auxil.∂θ∂t -= SOIL_BULK.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[i  ].auxil.δz;
        SOILS[i+1].auxil.∂θ∂t += SOIL_BULK.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[i+1].auxil.δz;
    end;

    return nothing
end;

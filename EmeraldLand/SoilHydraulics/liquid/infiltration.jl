# This file contains function related to the soil water infiltration

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-06: add function soil water infiltration
#
#######################################################################################################################################################################################################
"""

    soil_water_infiltration!(spac::BulkSPAC{FT}) where {FT}

Update the marginal increase of soil water content per layer, given
- `spac` `BulkSPAC` SPAC

"""
function soil_water_infiltration!(spac::BulkSPAC{FT}) where {FT}
    soil_bulk = spac.soil_bulk;
    soils = spac.soils;

    # update the soil bulk auxiliary variables related to soil water and enrgy transfer between adjacent layers
    N = length(soil_bulk.auxil.k);
    for i in 1:N
        soil_bulk.auxil.k[i]        = 1 / (2 / soils[i].auxil.k + 2 / soils[i+1].auxil.k);
        soil_bulk.auxil.δψ[i]       = soils[i].auxil.ψ - soils[i+1].auxil.ψ + ρg_MPa(FT) * (soils[i].auxil.z - soils[i+1].auxil.z);
        soil_bulk.auxil.q[i]        = soil_bulk.auxil.k[i] * soil_bulk.auxil.δψ[i];
        soil_bulk.auxil.λ_layers[i] = 1 / (2 / soils[i].auxil.λ_soil_water + 2 / soils[i+1].auxil.λ_soil_water);
        soil_bulk.auxil.δt[i]       = soils[i].auxil.t - soils[i+1].auxil.t;
        soil_bulk.auxil.q_layers[i] = soil_bulk.auxil.λ_layers[i] * soil_bulk.auxil.δt[i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (soil_bulk.auxil.q[i] > 0) && (soils[i+1].state.θ >= soils[i+1].state.vc.Θ_SAT)
            soil_bulk.auxil.q[i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (soil_bulk.auxil.q[i] < 0) && (soils[i].state.θ >= soils[i].state.vc.Θ_SAT)
            soil_bulk.auxil.q[i] = 0;
        end;

        # if both layers are oversaturated, move the oversaturated part from lower layer to upper layer
        if (soils[i].state.θ >= soils[i].state.vc.Θ_SAT) && (soils[i+1].state.θ > soils[i+1].state.vc.Θ_SAT)
            soil_bulk.auxil.q[i] = -1 * (soils[i+1].state.θ - soils[i+1].state.vc.Θ_SAT) * soils[i+1].auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);
        end;
    end;

    # update k, δψ, and flow rate among layers
    soils[1].auxil.∂θ∂t += spac.meteo.rain * M_H₂O(FT) / ρ_H₂O(FT) / soils[1].auxil.δz;

    N = length(soil_bulk.auxil.q);
    for i in 1:N
        soils[i  ].auxil.∂θ∂t -= soil_bulk.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / soils[i  ].auxil.δz;
        soils[i+1].auxil.∂θ∂t += soil_bulk.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / soils[i+1].auxil.δz;
    end;

    return nothing
end;

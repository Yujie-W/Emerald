# This file contains function related to the soil water infiltration

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-06: add function soil water infiltration
#     2024-Feb-29: add controller to make sure when both layers are very dry (potential < -50 MPa), there is no water exchange between them
#
#######################################################################################################################################################################################################
"""

    soil_water_infiltration!(spac::BulkSPAC{FT}) where {FT}

Update the marginal increase of soil water content per layer, given
- `spac` `BulkSPAC` SPAC

"""
function soil_water_infiltration!(spac::BulkSPAC{FT}) where {FT}
    sbulk = spac.soil_bulk;
    soils = spac.soils;

    # update the soil bulk auxiliary variables related to soil water and enrgy transfer between adjacent layers
    N = length(sbulk.auxil.k);
    for i in 1:N
        sbulk.auxil.k[i] = 1 / (2 / soils[i].s_aux.k + 2 / soils[i+1].s_aux.k);
        sbulk.auxil.δψ[i] = soils[i].s_aux.ψ - soils[i+1].s_aux.ψ + ρg_MPa(FT) * (soils[i].t_aux.z - soils[i+1].t_aux.z);
        sbulk.auxil.q[i] = sbulk.auxil.k[i] * sbulk.auxil.δψ[i];
        sbulk.auxil.λ_layers[i] = 1 / (2 / soils[i].s_aux.λ_soil_water + 2 / soils[i+1].s_aux.λ_soil_water);
        sbulk.auxil.δt[i] = soils[i].s_aux.t - soils[i+1].s_aux.t;
        sbulk.auxil.q_layers[i] = sbulk.auxil.λ_layers[i] * sbulk.auxil.δt[i];

        # if both layers are very dry, set the flow to 0
        if soils[i].s_aux.ψ < -50 && soils[i+1].s_aux.ψ < -50
            sbulk.auxil.δψ[i] = 0;
            sbulk.auxil.q[i] = 0;
        end;

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (sbulk.auxil.q[i] > 0) && (soils[i+1].state.θ >= soils[i+1].trait.vc.Θ_SAT)
            sbulk.auxil.q[i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (sbulk.auxil.q[i] < 0) && (soils[i].state.θ >= soils[i].trait.vc.Θ_SAT)
            sbulk.auxil.q[i] = 0;
        end;

        # if both layers are oversaturated, move the oversaturated part from lower layer to upper layer
        if (soils[i].state.θ >= soils[i].trait.vc.Θ_SAT) && (soils[i+1].state.θ > soils[i+1].trait.vc.Θ_SAT)
            sbulk.auxil.q[i] = -1 * (soils[i+1].state.θ - soils[i+1].trait.vc.Θ_SAT) * soils[i+1].t_aux.δz * ρ_H₂O(FT) / M_H₂O(FT);
        end;
    end;

    # update k, δψ, and flow rate among layers
    soils[1].auxil.∂θ∂t += spac.meteo.rain * M_H₂O(FT) / ρ_H₂O(FT) / soils[1].t_aux.δz;

    N = length(sbulk.auxil.q);
    for i in 1:N
        soils[i  ].auxil.∂θ∂t -= sbulk.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / soils[i  ].t_aux.δz;
        soils[i+1].auxil.∂θ∂t += sbulk.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / soils[i+1].t_aux.δz;
    end;

    return nothing
end;

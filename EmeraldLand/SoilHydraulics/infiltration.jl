#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-29: copy function out of soil_budget!
#     2023-Jun-30: add DEBUG mode to check for NaNs
#     2023-Jul-06: add info into DEBUG code block
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-07: add integrators for soil water budget
#     2023-Sep-11: rename ALLOW_SOIL_EVAPORATION to ENABLE_SOIL_EVAPORATION
#
#######################################################################################################################################################################################################
"""
#
    soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the infiltration rate between soil layers, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

#
    soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update soil water content and energy per layer, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
function soil_infiltration! end

soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; DEBUG, DIM_SOIL) = config;
    (; METEO, SOIL_BULK, SOILS) = spac;

    # update soil k, kd, ψ, and λ_thermal for each soil layer (0.5 for tortuosity factor)
    for soil in SOILS
        soil.auxil.k         = relative_hydraulic_conductance(soil.state.vc, soil.state.θ) * soil.state.vc.K_MAX * relative_viscosity(soil.auxil.t) / soil.auxil.δz;
        soil.auxil.ψ         = soil_ψ_25(soil.state.vc, soil.state.θ; oversaturation = true) * relative_surface_tension(soil.auxil.t);
        soil.auxil.kd        = 0.5 * max(0, soil.state.vc.Θ_SAT - soil.state.θ) / soil.auxil.δz;
        soil.auxil.λ_thermal = (soil.state.λ_thermal + soil.state.θ * Λ_THERMAL_H₂O(FT)) / soil.auxil.δz;
        soil.auxil.∂e∂t      = 0;
        soil.auxil.∂n∂t     .= 0;
        soil.auxil.∂θ∂t      = 0;
    end;

    # update k, δψ, and flow rate among layers
    SOILS[1].auxil.∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[1].auxil.δz;
    SOILS[1].auxil.∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    SOILS[1].auxil.∂e∂t += SOIL_BULK.auxil.r_net_lw + SOIL_BULK.auxil.r_net_sw;

    for i in 1:DIM_SOIL-1
        SOIL_BULK.auxil.k[i]         = 1 / (2 / SOILS[i].auxil.k + 2 / SOILS[i+1].auxil.k);
        SOIL_BULK.auxil.δψ[i]        = SOILS[i].auxil.ψ - SOILS[i+1].auxil.ψ + ρg_MPa(FT) * (SOILS[i].auxil.z - SOILS[i+1].auxil.z);
        SOIL_BULK.auxil.q[i]         = SOIL_BULK.auxil.k[i] * SOIL_BULK.auxil.δψ[i];
        SOIL_BULK.auxil.λ_thermal[i] = 1 / (2 / SOILS[i].auxil.λ_thermal + 2 / SOILS[i+1].auxil.λ_thermal);
        SOIL_BULK.auxil.δt[i]        = SOILS[i].auxil.t - SOILS[i+1].auxil.t;
        SOIL_BULK.auxil.q_thermal[i] = SOIL_BULK.auxil.λ_thermal[i] * SOIL_BULK.auxil.δt[i];

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

        SOILS[i  ].auxil.∂θ∂t -= SOIL_BULK.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[i].auxil.δz;
        SOILS[i+1].auxil.∂θ∂t += SOIL_BULK.auxil.q[i] * M_H₂O(FT) / ρ_H₂O(FT) / SOILS[i+1].auxil.δz;
        SOILS[i  ].auxil.∂e∂t -= SOIL_BULK.auxil.q_thermal[i];
        SOILS[i+1].auxil.∂e∂t += SOIL_BULK.auxil.q_thermal[i];
        SOILS[i  ].auxil.∂e∂t -= SOIL_BULK.auxil.q[i] * CP_L_MOL(FT) * SOILS[i].auxil.t;
        SOILS[i+1].auxil.∂e∂t += SOIL_BULK.auxil.q[i] * CP_L_MOL(FT) * SOILS[i].auxil.t;
    end;

    return nothing
);

soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; ENABLE_SOIL_EVAPORATION, DEBUG) = config;
    (; SOILS) = spac;

    # run the water transport (condensation + mass flow)
    for i in eachindex(SOILS)
        soil = SOILS[i];

        # account for evaporation and condensation to/from the air space
        if ENABLE_SOIL_EVAPORATION
            _ps = saturation_vapor_pressure(soil.t, soil.auxil.ψ * 1000000);
            _δθ_v = (soil.state.ns[3] / soil.auxil.δz - _ps * max(0, soil.state.vc.Θ_SAT - soil.θ) / (GAS_R(FT) * soil.t)) * M_H₂O(FT) / ρ_H₂O(FT);

            soil.state.θ += _δθ_v;
            soil.state.Σe += _δθ_v * ρ_H₂O(FT) * CP_V(FT) * soil.t; # this energy is transferred from/to air, so use CP_V
            soil.state.Σe += _δθ_v * ρ_H₂O(FT) * latent_heat_vapor(soil.t);
        end;

        # account for mass flow
        soil.state.θ += soil.auxil.∂θ∂t * δt;
        soil.state.Σe += soil.auxil.∂e∂t * δt / soil.auxil.δz;
    end;

    return nothing
);

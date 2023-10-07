# This file contains functions to calaculate the trace gas diffusion in the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-07: add function to compute the diffusion of trace gasses in the soil
#
#######################################################################################################################################################################################################
"""

    trace_gas_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the diffusion of trace gasses in the soil, given
- `config` `SPACConfiguration` type configuration
- `spac` `MultiLayerSPAC` type SPAC

"""
function trace_gas_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; AIR, SOIL_BULK, SOILS) = spac;
    (; DIM_SOIL, TRACE_AIR, TRACE_CH₄, TRACE_CO₂, TRACE_H₂O, TRACE_N₂, TRACE_O₂) = config;

    #
    # diffusion from top soil to the first air layer
    #     - diffusion path length is half of the thickness of the first soil layer
    #     - if there is no air in the first soil layer (saturated soil), there will be no diffusion for other trace gasses except for H₂O
    #     - to enable water to continue to evaporate, we add a small amount of volume in every soil layer (for H₂O diffusion only, note to update ns[3] accordingly)
    # TODO:
    #     - the impact from wind speed on the boundary layer thickness need to be accounted for
    #
    air = AIR[1];
    rt_air = GAS_R(FT) * air.t;
    rel_kd = 2 * SOILS[1].auxil.kd;
    rel_kv = 2 * SOILS[1].auxil.kv;
    v_gas = max(0, SOILS[1].state.vc.Θ_SAT - SOILS[1].state.θ) * SOILS[1].auxil.δz;
    if v_gas > 0
        SOIL_BULK.auxil.dndt[1,1] = rel_kd * diffusive_coefficient(SOILS[1].auxil.t, TRACE_CH₄, TRACE_AIR) * (SOILS[1].state.ns[1] / v_gas - air.p_CH₄ / rt_air)
        SOIL_BULK.auxil.dndt[1,2] = rel_kd * diffusive_coefficient(SOILS[1].auxil.t, TRACE_CO₂, TRACE_AIR) * (SOILS[1].state.ns[2] / v_gas - air.p_CO₂ / rt_air)
        SOIL_BULK.auxil.dndt[1,4] = rel_kd * diffusive_coefficient(SOILS[1].auxil.t, TRACE_N₂ , TRACE_AIR) * (SOILS[1].state.ns[4] / v_gas - air.p_N₂  / rt_air)
        SOIL_BULK.auxil.dndt[1,5] = rel_kd * diffusive_coefficient(SOILS[1].auxil.t, TRACE_O₂ , TRACE_AIR) * (SOILS[1].state.ns[5] / v_gas - air.p_O₂  / rt_air)
    end;
    v_vap = v_gas + FT(0.01);
    SOIL_BULK.auxil.dndt[1,3] = rel_kv * diffusive_coefficient(SOILS[1].auxil.t, TRACE_H₂O, TRACE_AIR) * (SOILS[1].state.ns[3] / v_vap - air.p_H₂O / rt_air)

    #
    # Diffusion among soi layers
    #     - diffusion path length is the sum of the half thickness of the two layers
    #     - mean diffusion rate is the mean of the two layers 1/k = 1/k1 + 1/k2
    #
    for i in 1:DIM_SOIL-1
        # for gas 1,2,4,5
        v_gas_i = max(0, SOILS[i  ].state.vc.Θ_SAT - SOILS[i  ].state.θ) * SOILS[i  ].auxil.δz;
        v_gas_j = max(0, SOILS[i+1].state.vc.Θ_SAT - SOILS[i+1].state.θ) * SOILS[i+1].auxil.δz;
        if v_gas_i > 0 && v_gas_j > 0
            kdi = diffusive_coefficient(SOILS[i  ].auxil.t, TRACE_CH₄, TRACE_AIR) * SOILS[i  ].auxil.kd;
            kdj = diffusive_coefficient(SOILS[i+1].auxil.t, TRACE_CH₄, TRACE_AIR) * SOILS[i+1].auxil.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (SOILS[i].state.ns[1] / v_gas_i - SOILS[i+1].state.ns[1] / v_gas_j);
            SOIL_BULK.auxil.dndt[i+1,1] = -drate;

            kdi = diffusive_coefficient(SOILS[i  ].auxil.t, TRACE_CO₂, TRACE_AIR) * SOILS[i  ].auxil.kd;
            kdj = diffusive_coefficient(SOILS[i+1].auxil.t, TRACE_CO₂, TRACE_AIR) * SOILS[i+1].auxil.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (SOILS[i].state.ns[2] / v_gas_i - SOILS[i+1].state.ns[2] / v_gas_j);
            SOIL_BULK.auxil.dndt[i+1,2] = -drate;

            kdi = diffusive_coefficient(SOILS[i  ].auxil.t, TRACE_N₂, TRACE_AIR) * SOILS[i  ].auxil.kd;
            kdj = diffusive_coefficient(SOILS[i+1].auxil.t, TRACE_N₂, TRACE_AIR) * SOILS[i+1].auxil.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (SOILS[i].state.ns[4] / v_gas_i - SOILS[i+1].state.ns[4] / v_gas_j);
            SOIL_BULK.auxil.dndt[i+1,4] = -drate;

            kdi = diffusive_coefficient(SOILS[i  ].auxil.t, TRACE_O₂, TRACE_AIR) * SOILS[i  ].auxil.kd;
            kdj = diffusive_coefficient(SOILS[i+1].auxil.t, TRACE_O₂, TRACE_AIR) * SOILS[i+1].auxil.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (SOILS[i].state.ns[5] / v_gas_i - SOILS[i+1].state.ns[5] / v_gas_j);
            SOIL_BULK.auxil.dndt[i+1,5] = -drate;
        end;

        # for water vapor
        v_vap_i = v_gas_i + FT(0.01);
        v_vap_j = v_gas_j + FT(0.01);
        kdi = diffusive_coefficient(SOILS[i  ].auxil.t, TRACE_H₂O, TRACE_AIR) * SOILS[i  ].auxil.kv;
        kdj = diffusive_coefficient(SOILS[i+1].auxil.t, TRACE_H₂O, TRACE_AIR) * SOILS[i+1].auxil.kv;
        rel_kv = 2 * kdi * kdj / (kdi + kdj);
        drate = rel_kv * (SOILS[i].state.ns[3] / v_vap_i - SOILS[i+1].state.ns[3] / v_vap_j);
        SOIL_BULK.auxil.dndt[i+1,3] = -drate;
    end;

    # update the diffusion from top soil to the first air layer
    for j in 1:5
        SOILS[1].auxil.∂n∂t[j] -= SOIL_BULK.auxil.dndt[1,j];
    end;

    # update the diffusion among soil layers
    for i in 1:DIM_SOIL-1
        for j in 1:5
            SOILS[i  ].auxil.∂n∂t[j] += SOIL_BULK.auxil.dndt[i+1,j];
            SOILS[i+1].auxil.∂n∂t[j] -= SOIL_BULK.auxil.dndt[i+1,j];
        end;
    end;

    return nothing
end;

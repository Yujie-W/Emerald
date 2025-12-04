# This file contains functions to calaculate the trace gas diffusion in the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-07: add function to compute the diffusion of trace gasses in the soil
#     2025-Jun-05: account for ice volume in the gas diffusion
#
#######################################################################################################################################################################################################
"""

    trace_gas_diffusion!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Compute the diffusion of trace gasses in the soil, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` type SPAC

"""
function trace_gas_diffusion!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    (; TRACE_AIR, TRACE_CH₄, TRACE_CO₂, TRACE_H₂O, TRACE_N₂, TRACE_O₂) = config;
    sbulk = spac.soil_bulk;
    soils = spac.soils;
    n_soil = length(soils);

    #
    # diffusion from top soil to the first air layer
    #     - diffusion path length is half of the thickness of the first soil layer
    #     - if there is no air in the first soil layer (saturated soil), there will be no diffusion for other trace gasses except for H₂O
    #     - to enable water to continue to evaporate, we add a small amount of volume in every soil layer (for H₂O diffusion only, note to update ns[3] accordingly)
    # TODO:
    #     - the impact from wind speed on the boundary layer thickness need to be accounted for
    #     - the impact of ice volume change not yet accounted for
    #
    air = spac.airs[1];
    rt_air = GAS_R(FT) * air.s_aux.t;
    rel_kd = 2 * soils[1].s_aux.kd;
    rel_kv = 2 * soils[1].s_aux.kv;
    v_gas = max(0, soils[1].trait.vc.Θ_SAT - soils[1].state.θ - soils[1].state.θ_ice) * soils[1].t_aux.δz;
    if v_gas > 0
        sbulk.auxil.dndt[1,1] = rel_kd * diffusive_coefficient(soils[1].s_aux.t, TRACE_CH₄, TRACE_AIR) * (soils[1].state.ns[1] / v_gas - air.s_aux.ps[1] / rt_air);
        sbulk.auxil.dndt[1,2] = rel_kd * diffusive_coefficient(soils[1].s_aux.t, TRACE_CO₂, TRACE_AIR) * (soils[1].state.ns[2] / v_gas - air.s_aux.ps[2] / rt_air);
        sbulk.auxil.dndt[1,4] = rel_kd * diffusive_coefficient(soils[1].s_aux.t, TRACE_N₂ , TRACE_AIR) * (soils[1].state.ns[4] / v_gas - air.s_aux.ps[4] / rt_air);
        sbulk.auxil.dndt[1,5] = rel_kd * diffusive_coefficient(soils[1].s_aux.t, TRACE_O₂ , TRACE_AIR) * (soils[1].state.ns[5] / v_gas - air.s_aux.ps[5] / rt_air);
    end;
    v_vap = v_gas + FT(0.01);
    sbulk.auxil.dndt[1,3] = rel_kv * diffusive_coefficient(soils[1].s_aux.t, TRACE_H₂O, TRACE_AIR) * (soils[1].state.ns[3] / v_vap - air.s_aux.ps[3] / rt_air);

    #
    # Diffusion among soi layers
    #     - diffusion path length is the sum of the half thickness of the two layers
    #     - mean diffusion rate is the mean of the two layers 1/k = 1/k1 + 1/k2
    #
    for i in 1:n_soil-1
        # for gas 1,2,4,5
        v_gas_i = max(0, soils[i  ].trait.vc.Θ_SAT - soils[i  ].state.θ - soils[i  ].state.θ_ice) * soils[i  ].t_aux.δz;
        v_gas_j = max(0, soils[i+1].trait.vc.Θ_SAT - soils[i+1].state.θ - soils[i+1].state.θ_ice) * soils[i+1].t_aux.δz;
        if v_gas_i > 0 && v_gas_j > 0
            kdi = diffusive_coefficient(soils[i  ].s_aux.t, TRACE_CH₄, TRACE_AIR) * soils[i  ].s_aux.kd;
            kdj = diffusive_coefficient(soils[i+1].s_aux.t, TRACE_CH₄, TRACE_AIR) * soils[i+1].s_aux.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (soils[i].state.ns[1] / v_gas_i - soils[i+1].state.ns[1] / v_gas_j);
            sbulk.auxil.dndt[i+1,1] = -drate;

            kdi = diffusive_coefficient(soils[i  ].s_aux.t, TRACE_CO₂, TRACE_AIR) * soils[i  ].s_aux.kd;
            kdj = diffusive_coefficient(soils[i+1].s_aux.t, TRACE_CO₂, TRACE_AIR) * soils[i+1].s_aux.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (soils[i].state.ns[2] / v_gas_i - soils[i+1].state.ns[2] / v_gas_j);
            sbulk.auxil.dndt[i+1,2] = -drate;

            kdi = diffusive_coefficient(soils[i  ].s_aux.t, TRACE_N₂, TRACE_AIR) * soils[i  ].s_aux.kd;
            kdj = diffusive_coefficient(soils[i+1].s_aux.t, TRACE_N₂, TRACE_AIR) * soils[i+1].s_aux.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (soils[i].state.ns[4] / v_gas_i - soils[i+1].state.ns[4] / v_gas_j);
            sbulk.auxil.dndt[i+1,4] = -drate;

            kdi = diffusive_coefficient(soils[i  ].s_aux.t, TRACE_O₂, TRACE_AIR) * soils[i  ].s_aux.kd;
            kdj = diffusive_coefficient(soils[i+1].s_aux.t, TRACE_O₂, TRACE_AIR) * soils[i+1].s_aux.kd;
            rel_kd = 2 * kdi * kdj / (kdi + kdj);
            drate = rel_kd * (soils[i].state.ns[5] / v_gas_i - soils[i+1].state.ns[5] / v_gas_j);
            sbulk.auxil.dndt[i+1,5] = -drate;
        end;

        # for water vapor
        v_vap_i = v_gas_i + FT(0.01);
        v_vap_j = v_gas_j + FT(0.01);
        kdi = diffusive_coefficient(soils[i  ].s_aux.t, TRACE_H₂O, TRACE_AIR) * soils[i  ].s_aux.kv;
        kdj = diffusive_coefficient(soils[i+1].s_aux.t, TRACE_H₂O, TRACE_AIR) * soils[i+1].s_aux.kv;
        rel_kv = 2 * kdi * kdj / (kdi + kdj);
        drate = rel_kv * (soils[i].state.ns[3] / v_vap_i - soils[i+1].state.ns[3] / v_vap_j);
        sbulk.auxil.dndt[i+1,3] = -drate;
    end;

    # update the diffusion from top soil to the first air layer
    for j in 1:5
        soils[1].auxil.∂n∂t[j] -= sbulk.auxil.dndt[1,j];
    end;

    # update the diffusion among soil layers
    for i in 1:n_soil-1
        for j in 1:5
            soils[i  ].auxil.∂n∂t[j] += sbulk.auxil.dndt[i+1,j];
            soils[i+1].auxil.∂n∂t[j] -= sbulk.auxil.dndt[i+1,j];
        end;
    end;

    return nothing
end;

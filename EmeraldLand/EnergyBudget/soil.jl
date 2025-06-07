# This file contains function to calculate energy budgets of the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-06: add function soil energy flow
#     2023-Oct-07: add calculations for diffusion related energy transfer
#     2025-Jun-05: all soil enery flow calculations are relative to the triple point temperature (T₀)
#     2026-Jun-06: account for sensible heat transfer between soil and air (assumed 0.2 mol m⁻² s⁻¹ conductance)
# Todos
#     TODO: compute the sensible heat transfer between soil and air using wind speed of the first air layer
#
#######################################################################################################################################################################################################
"""

    soil_energy_flow!(spac::BulkSPAC{FT}) where {FT}

Update the marginal increase of soil energy per layer, given
- `spac` `BulkSPAC` SPAC

Note that energy associated with water movement between soil and root is accounted for in the root energy budget function.

"""
function soil_energy_flow!(spac::BulkSPAC{FT}) where {FT}
    air = spac.airs[1];
    meteo = spac.meteo;
    sbulk = spac.soil_bulk;
    soils = spac.soils;

    # add rain and radiation energy into the first layer
    soils[1].auxil.∂e∂t += meteo.rain * CP_L_MOL(FT) * (meteo.t_precip - T₀(FT));
    soils[1].auxil.∂e∂t += sbulk.auxil.r_net_lw + sbulk.auxil.r_net_sw;

    # water and energy exchange among layers
    N = length(sbulk.auxil.q);
    for i in 1:N
        soils[i  ].auxil.∂e∂t -= sbulk.auxil.q_layers[i];
        soils[i+1].auxil.∂e∂t += sbulk.auxil.q_layers[i];

        # if flow from upper to lower is positive, use temperature from upper layer; otherwise, use temperature from lower layer
        t = sbulk.auxil.q[i] >= 0 ? soils[i].s_aux.t : soils[i+1].s_aux.t;
        soils[i  ].auxil.∂e∂t -= sbulk.auxil.q[i] * CP_L_MOL(FT) * (t - T₀(FT));
        soils[i+1].auxil.∂e∂t += sbulk.auxil.q[i] * CP_L_MOL(FT) * (t - T₀(FT));
    end;

    # energy transfer related to gas diffusion (CH₄, CO₂, H₂O, N₂, O₂)
    # update the energy for diffusion from top soil to the first air layer
    δn1 = sbulk.auxil.dndt[1,3];
    δn4 = sbulk.auxil.dndt[1,1] + sbulk.auxil.dndt[1,2] + sbulk.auxil.dndt[1,4] + sbulk.auxil.dndt[1,5];
    t = δn1 > 0 ? soils[1].s_aux.t : air.s_aux.t;
    soils[1].auxil.∂e∂t -= δn1 * CP_V_MOL(FT) * (t - T₀(FT));
    t = δn4 > 0 ? soils[1].s_aux.t : air.s_aux.t;
    soils[1].auxil.∂e∂t -= δn4 * CP_D_MOL(FT) * (t - T₀(FT));

    # update the diffusion among soil layers
    N = length(sbulk.auxil.q);
    for i in 1:N
        δn1 = sbulk.auxil.dndt[i+1,3];
        δn4 = sbulk.auxil.dndt[i+1,1] + sbulk.auxil.dndt[i+1,2] + sbulk.auxil.dndt[i+1,4] + sbulk.auxil.dndt[i+1,5];
        t = δn1 > 0 ? soils[i+1].s_aux.t : soils[i].s_aux.t;
        soils[i  ].auxil.∂e∂t += δn1 * CP_V_MOL(FT) * (t - T₀(FT));
        soils[i+1].auxil.∂e∂t -= δn1 * CP_V_MOL(FT) * (t - T₀(FT));
        t = δn4 > 0 ? soils[i+1].s_aux.t : soils[i].s_aux.t;
        soils[i  ].auxil.∂e∂t += δn4 * CP_D_MOL(FT) * (t - T₀(FT));
        soils[i+1].auxil.∂e∂t -= δn4 * CP_D_MOL(FT) * (t - T₀(FT));
    end;

    # TODO: assume a 0.2 mol m⁻² s⁻¹ heat conductance between the first air layer and the first soil layer
    soils[1].auxil.∂e∂t -= 0.2 * CP_D_MOL(FT) * (soils[1].s_aux.t - air.s_aux.t);

    return nothing
end;

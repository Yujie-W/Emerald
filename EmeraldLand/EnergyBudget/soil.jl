# This file contains function to calculate energy budgets of the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-06: add function soil energy flow
#     2023-Oct-07: add calculations for diffusion related energy transfer
#
#######################################################################################################################################################################################################
"""

    soil_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}

Update the marginal increase of soil energy per layer, given
- `spac` `MultiLayerSPAC` SPAC

Note that energy associated with water movement between soil and root is accounted for in the root energy budget function.

"""
function soil_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}
    (; AIR, METEO, SOIL_BULK, SOILS) = spac;

    # add rain and radiation energy into the first layer
    SOILS[1].auxil.∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    SOILS[1].auxil.∂e∂t += SOIL_BULK.auxil.r_net_lw + SOIL_BULK.auxil.r_net_sw;

    # water and energy exchange among layers
    N = length(SOIL_BULK.auxil.q);
    for i in 1:N
        SOILS[i  ].auxil.∂e∂t -= SOIL_BULK.auxil.q_layers[i];
        SOILS[i+1].auxil.∂e∂t += SOIL_BULK.auxil.q_layers[i];

        # if flow from upper to lower is positive, use temperature from upper layer; otherwise, use temperature from lower layer
        t = SOIL_BULK.auxil.q[i] >= 0 ? SOILS[i].auxil.t : SOILS[i+1].auxil.t;
        SOILS[i  ].auxil.∂e∂t -= SOIL_BULK.auxil.q[i] * CP_L_MOL(FT) * t;
        SOILS[i+1].auxil.∂e∂t += SOIL_BULK.auxil.q[i] * CP_L_MOL(FT) * t;
    end;

    # energy transfer related to gas diffusion (CH₄, CO₂, H₂O, N₂, O₂)
    # update the energy for diffusion from top soil to the first air layer
    δn1 = SOIL_BULK.auxil.dndt[1,3];
    δn4 = SOIL_BULK.auxil.dndt[1,1] + SOIL_BULK.auxil.dndt[1,2] + SOIL_BULK.auxil.dndt[1,4] + SOIL_BULK.auxil.dndt[1,5];
    t = δn1 > 0 ? SOILS[1].auxil.t : AIR[1].t;
    SOILS[1].auxil.∂e∂t -= δn1 * CP_V_MOL(FT) * t;
    t = δn4 > 0 ? SOILS[1].auxil.t : AIR[1].t;
    SOILS[1].auxil.∂e∂t -= δn4 * CP_D_MOL(FT) * t;

    # update the diffusion among soil layers
    N = length(SOIL_BULK.auxil.q);
    for i in 1:N
        δn1 = SOIL_BULK.auxil.dndt[i+1,3];
        δn4 = SOIL_BULK.auxil.dndt[i+1,1] + SOIL_BULK.auxil.dndt[i+1,2] + SOIL_BULK.auxil.dndt[i+1,4] + SOIL_BULK.auxil.dndt[i+1,5];
        t = δn1 > 0 ? SOILS[i+1].auxil.t : SOILS[i].auxil.t;
        SOILS[i  ].auxil.∂e∂t += δn1 * CP_V_MOL(FT) * t;
        SOILS[i+1].auxil.∂e∂t -= δn1 * CP_V_MOL(FT) * t;
        t = δn4 > 0 ? SOILS[i+1].auxil.t : SOILS[i].auxil.t;
        SOILS[i  ].auxil.∂e∂t += δn4 * CP_D_MOL(FT) * t;
        SOILS[i+1].auxil.∂e∂t -= δn4 * CP_D_MOL(FT) * t;
    end;

    return nothing
end;

# This file contains function to calculate energy budgets of the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-06: add function soil energy flow
#
#######################################################################################################################################################################################################
"""

    soil_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}

Update the marginal increase of soil energy per layer, given
- `spac` `MultiLayerSPAC` SPAC

"""
function soil_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}
    (; METEO, SOIL_BULK, SOILS) = spac;

    # add rain and radiation energy into the first layer
    SOILS[1].auxil.∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    SOILS[1].auxil.∂e∂t += SOIL_BULK.auxil.r_net_lw + SOIL_BULK.auxil.r_net_sw;

    N = length(SOIL_BULK.auxil.q);
    for i in 1:N
        SOILS[i  ].auxil.∂e∂t -= SOIL_BULK.auxil.q_layers[i];
        SOILS[i+1].auxil.∂e∂t += SOIL_BULK.auxil.q_layers[i];

        # if flow from upper to lower is positive, use temperature from upper layer; otherwise, use temperature from lower layer
        t = SOIL_BULK.auxil.q[i] >= 0 ? SOILS[i].auxil.t : SOILS[i+1].auxil.t;
        SOILS[i  ].auxil.∂e∂t -= SOIL_BULK.auxil.q[i] * CP_L_MOL(FT) * t;
        SOILS[i+1].auxil.∂e∂t += SOIL_BULK.auxil.q[i] * CP_L_MOL(FT) * t;
    end;

    return nothing
end;

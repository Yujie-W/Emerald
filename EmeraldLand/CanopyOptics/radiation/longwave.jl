# This file contains function to run canopy longwave radiation simulations

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-11: add function longwave_radiation!
#     2023-Oct-14: if LAI <= 0, run soil longwave radiation only
#
#######################################################################################################################################################################################################
"""

    longwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Run longwave radiation simulations, given
- `config` SPAC configuration
- `spac` SPAC

"""
function longwave_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    (; DIM_LAYER) = config;
    can_struct = spac.canopy.structure;
    leaves = spac.plant.leaves;
    meteo = spac.meteo;
    soil_bulk = spac.soil_bulk;
    top_soil = spac.soils[1];

    if can_struct.state.lai <= 0 && can_struct.state.sai <= 0
        # 1. compute longwave radiation out from the leaves and soil
        can_struct.auxil.lw_layer .= 0;
        r_lw_soil = K_STEFAN(FT) * (1 - soil_bulk.auxil.ρ_lw) * top_soil.auxil.t ^ 4;

        # 2. account for the longwave emission from bottom to up
        can_struct.auxil.emitꜜ .= 0;
        can_struct.auxil.emitꜛ .= r_lw_soil;

        # 3. account for the longwave emission from up to bottom
        can_struct.auxil.lwꜜ .= meteo.rad_lw;
        can_struct.auxil.lwꜛ .= meteo.rad_lw * soil_bulk.auxil.ρ_lw + r_lw_soil;

        # 4. compute the net longwave radiation per canopy layer and soil
        can_struct.auxil.r_net_lw .= 0;
        soil_bulk.auxil.r_net_lw = meteo.rad_lw * (1 - soil_bulk.auxil.ρ_lw) - r_lw_soil;

        return nothing
    end;

    # 1. compute longwave radiation out from the leaves and soil
    for i in 1:DIM_LAYER
        leaf = leaves[DIM_LAYER + 1 - i];
        can_struct.auxil.lw_layer[i] = K_STEFAN(FT) * can_struct.auxil.ϵ_lw_layer[i] * leaf.energy.auxil.t ^ 4;
    end;
    r_lw_soil = K_STEFAN(FT) * (1 - soil_bulk.auxil.ρ_lw) * top_soil.auxil.t ^ 4;

    # 2. account for the longwave emission from bottom to up
    can_struct.auxil.emitꜛ[end] = r_lw_soil;
    for i in DIM_LAYER:-1:1
        r = can_struct.auxil.ρ_lw_layer[i];
        t = can_struct.auxil.τ_lw_layer[i];
        r_j = can_struct.auxil.ρ_lw[i+1];

        denom = 1 - r * r_j;
        can_struct.auxil.emitꜜ[i] = (can_struct.auxil.emitꜛ[i+1] * r + can_struct.auxil.lw_layer[i]) / denom;
        can_struct.auxil.emitꜛ[i] = (can_struct.auxil.emitꜜ[i] * r_j * t + can_struct.auxil.emitꜛ[i+1] * t) / denom + can_struct.auxil.lw_layer[i];
    end;

    # 3. account for the longwave emission from up to bottom
    can_struct.auxil.lwꜜ[1] = meteo.rad_lw;
    for i in 1:DIM_LAYER
        r_i = can_struct.auxil.ρ_lw[i];
        t_i = can_struct.auxil.τ_lw[i];

        can_struct.auxil.lwꜛ[i] = can_struct.auxil.lwꜜ[i] * r_i + can_struct.auxil.emitꜛ[i];
        can_struct.auxil.lwꜜ[i+1] = can_struct.auxil.lwꜜ[i] * t_i + can_struct.auxil.emitꜜ[i];
    end;
    can_struct.auxil.lwꜛ[end] = can_struct.auxil.lwꜜ[end] * soil_bulk.auxil.ρ_lw + r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    # TODO: stem energy balance
    for i in 1:DIM_LAYER
        can_struct.auxil.r_net_lw[i] = (can_struct.auxil.lwꜜ[i] + can_struct.auxil.lwꜛ[i+1]) * can_struct.auxil.ϵ_lw_layer[i] - 2 * can_struct.auxil.lw_layer[i];
        can_struct.auxil.r_net_lw[i] /= can_struct.state.δlai[i];
    end;

    soil_bulk.auxil.r_net_lw = can_struct.auxil.lwꜜ[end] * (1 - soil_bulk.auxil.ρ_lw) - r_lw_soil;

    return nothing
end;

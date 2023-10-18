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
    can_str = spac.canopy.structure;
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;
    meteo = spac.meteo;
    sbulk = spac.soil_bulk;
    top_soil = spac.soils[1];

    if can_str.state.lai <= 0 && can_str.state.sai <= 0
        # 1. compute longwave radiation out from the leaves and soil
        can_str.auxil.lw_layer .= 0;
        r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * top_soil.auxil.t ^ 4;

        # 2. account for the longwave emission from bottom to up
        can_str.auxil.emitꜜ .= 0;
        can_str.auxil.emitꜛ .= r_lw_soil;

        # 3. account for the longwave emission from up to bottom
        can_str.auxil.lwꜜ .= meteo.rad_lw;
        can_str.auxil.lwꜛ .= meteo.rad_lw * sbulk.auxil.ρ_lw + r_lw_soil;

        # 4. compute the net longwave radiation per canopy layer and soil
        can_str.auxil.r_net_lw_leaf .= 0;
        can_str.auxil.r_net_lw_stem .= 0;
        sbulk.auxil.r_net_lw = meteo.rad_lw * (1 - sbulk.auxil.ρ_lw) - r_lw_soil;

        return nothing
    end;

    # 1. compute longwave radiation out from the leaves and soil
    for i in 1:DIM_LAYER
        leaf = leaves[DIM_LAYER + 1 - i];
        stem = branches[DIM_LAYER + 1 - i];
        can_str.auxil.lw_layer[i] = K_STEFAN(FT) * can_str.auxil.ϵ_lw_layer[i] * leaf.energy.auxil.t ^ 4;
    end;
    r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * top_soil.auxil.t ^ 4;

    # 2. account for the longwave emission from bottom to up
    can_str.auxil.emitꜛ[end] = r_lw_soil;
    for i in DIM_LAYER:-1:1
        r = can_str.auxil.ρ_lw_layer[i];
        t = can_str.auxil.τ_lw_layer[i];
        r_j = can_str.auxil.ρ_lw[i+1];

        denom = 1 - r * r_j;
        can_str.auxil.emitꜜ[i] = (can_str.auxil.emitꜛ[i+1] * r + can_str.auxil.lw_layer[i]) / denom;
        can_str.auxil.emitꜛ[i] = (can_str.auxil.emitꜜ[i] * r_j * t + can_str.auxil.emitꜛ[i+1] * t) / denom + can_str.auxil.lw_layer[i];
    end;

    # 3. account for the longwave emission from up to bottom
    can_str.auxil.lwꜜ[1] = meteo.rad_lw;
    for i in 1:DIM_LAYER
        r_i = can_str.auxil.ρ_lw[i];
        t_i = can_str.auxil.τ_lw[i];

        can_str.auxil.lwꜛ[i] = can_str.auxil.lwꜜ[i] * r_i + can_str.auxil.emitꜛ[i];
        can_str.auxil.lwꜜ[i+1] = can_str.auxil.lwꜜ[i] * t_i + can_str.auxil.emitꜜ[i];
    end;
    can_str.auxil.lwꜛ[end] = can_str.auxil.lwꜜ[end] * sbulk.auxil.ρ_lw + r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    # TODO: stem energy balance
    for i in 1:DIM_LAYER
        can_str.auxil.r_net_lw[i] = (can_str.auxil.lwꜜ[i] + can_str.auxil.lwꜛ[i+1]) * can_str.auxil.ϵ_lw_layer[i] - 2 * can_str.auxil.lw_layer[i];
        can_str.auxil.r_net_lw[i] /= can_str.state.δlai[i];
    end;

    sbulk.auxil.r_net_lw = can_str.auxil.lwꜜ[end] * (1 - sbulk.auxil.ρ_lw) - r_lw_soil;

    return nothing
end;

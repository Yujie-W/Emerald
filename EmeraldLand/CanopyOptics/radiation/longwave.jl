# This file contains function to run canopy longwave radiation simulations

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-11: add function longwave_radiation!
#     2023-Oct-14: if LAI <= 0, run soil longwave radiation only
#     2023-Oct-18: partition the energy between leaf and stem
#
#######################################################################################################################################################################################################
"""

    longwave_radiation!(spac::BulkSPAC{FT}) where {FT}

Run longwave radiation simulations, given
- `spac` SPAC

"""
function longwave_radiation!(spac::BulkSPAC{FT}) where {FT}
    branches = spac.plant.branches;
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    meteo = spac.meteo;
    sbulk = spac.soil_bulk;
    top_soil = spac.soils[1];
    n_layer = length(leaves);

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
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        stem = branches[ilf];
        # can_str.auxil.lw_layer[i] = K_STEFAN(FT) * can_str.auxil.ϵ_lw_layer[i] * leaf.energy.auxil.t ^ 4;
        f_leaf = can_str.state.δlai[irt] / (can_str.state.δlai[irt] + can_str.state.δsai[irt]);
        f_stem = 1 - f_leaf;
        can_str.auxil.lw_layer_leaf[irt] = leaf.energy.auxil.t ^ 4 * f_leaf * K_STEFAN(FT) * can_str.auxil.ϵ_lw_layer[irt];
        can_str.auxil.lw_layer_stem[irt] = stem.energy.auxil.t ^ 4 * f_stem * K_STEFAN(FT) * can_str.auxil.ϵ_lw_layer[irt];
        can_str.auxil.lw_layer[irt] = can_str.auxil.lw_layer_leaf[irt] + can_str.auxil.lw_layer_stem[irt];
    end;
    r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * top_soil.auxil.t ^ 4;

    # 2. account for the longwave emission from bottom to up
    can_str.auxil.emitꜛ[end] = r_lw_soil;
    for i in n_layer:-1:1
        r = can_str.auxil.ρ_lw_layer[i];
        t = can_str.auxil.τ_lw_layer[i];
        r_j = can_str.auxil.ρ_lw[i+1];

        denom = 1 - r * r_j;
        can_str.auxil.emitꜜ[i] = (can_str.auxil.emitꜛ[i+1] * r + can_str.auxil.lw_layer[i]) / denom;
        can_str.auxil.emitꜛ[i] = (can_str.auxil.emitꜜ[i] * r_j * t + can_str.auxil.emitꜛ[i+1] * t) / denom + can_str.auxil.lw_layer[i];
    end;

    # 3. account for the longwave emission from up to bottom
    can_str.auxil.lwꜜ[1] = meteo.rad_lw;
    for i in 1:n_layer
        r_i = can_str.auxil.ρ_lw[i];
        t_i = can_str.auxil.τ_lw[i];

        can_str.auxil.lwꜛ[i] = can_str.auxil.lwꜜ[i] * r_i + can_str.auxil.emitꜛ[i];
        can_str.auxil.lwꜜ[i+1] = can_str.auxil.lwꜜ[i] * t_i + can_str.auxil.emitꜜ[i];
    end;
    can_str.auxil.lwꜛ[end] = can_str.auxil.lwꜜ[end] * sbulk.auxil.ρ_lw + r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    for i in 1:n_layer
        # can_str.auxil.r_net_lw[i] = (can_str.auxil.lwꜜ[i] + can_str.auxil.lwꜛ[i+1]) * can_str.auxil.ϵ_lw_layer[i] - 2 * can_str.auxil.lw_layer[i];
        # can_str.auxil.r_net_lw[i] /= can_str.state.δlai[i];
        f_leaf = can_str.state.δlai[i] / (can_str.state.δlai[i] + can_str.state.δsai[i]);
        f_stem = 1 - f_leaf;
        can_str.auxil.r_net_lw_leaf[i] = (can_str.auxil.lwꜜ[i] + can_str.auxil.lwꜛ[i+1]) * can_str.auxil.ϵ_lw_layer[i] * f_leaf - 2 * can_str.auxil.lw_layer_leaf[i];
        can_str.auxil.r_net_lw_stem[i] = (can_str.auxil.lwꜜ[i] + can_str.auxil.lwꜛ[i+1]) * can_str.auxil.ϵ_lw_layer[i] * f_stem - 2 * can_str.auxil.lw_layer_stem[i];
    end;

    sbulk.auxil.r_net_lw = can_str.auxil.lwꜜ[end] * (1 - sbulk.auxil.ρ_lw) - r_lw_soil;

    return nothing
end;

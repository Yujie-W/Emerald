# This file contains function to run canopy longwave radiation simulations

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-11: add function longwave_radiation!
#
#######################################################################################################################################################################################################
"""

    longwave_radiation!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Run longwave radiation simulations, given
- `config` SPAC configuration
- `spac` SPAC

"""
function longwave_radiation!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; DIM_LAYER) = config;
    (; CANOPY, LEAVES, METEO, SOIL_BULK, SOILS) = spac;

    # 1. compute longwave radiation out from the leaves and soil
    for i in eachindex(LEAVES)
        CANOPY.structure.auxil.lw_layer[i] = K_STEFAN(FT) * CANOPY.structure.auxil.ϵ_lw_layer[i] * LEAVES[DIM_LAYER+1-i].energy.auxil.t ^ 4;
    end;
    r_lw_soil = K_STEFAN(FT) * (1 - SOIL_BULK.auxil.ρ_lw) * SOILS[1].auxil.t ^ 4;

    # 2. account for the longwave emission from bottom to up
    CANOPY.structure.auxil.emitꜛ[end] = r_lw_soil;
    for i in DIM_LAYER:-1:1
        r = CANOPY.structure.auxil.ρ_lw_layer[i];
        t = CANOPY.structure.auxil.τ_lw_layer[i];
        r_j = CANOPY.structure.auxil.ρ_lw[i+1];

        denom = 1 - r * r_j;
        CANOPY.structure.auxil.emitꜜ[i] = (CANOPY.structure.auxil.emitꜛ[i+1] * r + CANOPY.structure.auxil.lw_layer[i]) / denom;
        CANOPY.structure.auxil.emitꜛ[i] = (CANOPY.structure.auxil.emitꜜ[i] * r_j * t + CANOPY.structure.auxil.emitꜛ[i+1] * t) / denom + CANOPY.structure.auxil.lw_layer[i];
    end;

    # 3. account for the longwave emission from up to bottom
    CANOPY.structure.auxil.lwꜜ[1] = METEO.rad_lw;
    for i in 1:DIM_LAYER
        r_i = CANOPY.structure.auxil.ρ_lw[i];
        t_i = CANOPY.structure.auxil.τ_lw[i];

        CANOPY.structure.auxil.lwꜛ[i] = CANOPY.structure.auxil.lwꜜ[i] * r_i + CANOPY.structure.auxil.emitꜛ[i];
        CANOPY.structure.auxil.lwꜜ[i+1] = CANOPY.structure.auxil.lwꜜ[i] * t_i + CANOPY.structure.auxil.emitꜜ[i];
    end;
    CANOPY.structure.auxil.lwꜛ[end] = CANOPY.structure.auxil.lwꜜ[end] * SOIL_BULK.auxil.ρ_lw + r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    for i in 1:DIM_LAYER
        CANOPY.structure.auxil.r_net_lw[i] = (CANOPY.structure.auxil.lwꜜ[i] + CANOPY.structure.auxil.lwꜛ[i+1]) * CANOPY.structure.auxil.ϵ_lw_layer[i] - 2 * CANOPY.structure.auxil.lw_layer[i];
        CANOPY.structure.auxil.r_net_lw[i] /= CANOPY.structure.state.δlai[i];
    end;

    SOIL_BULK.auxil.r_net_lw = CANOPY.structure.auxil.lwꜜ[end] * (1 - SOIL_BULK.auxil.ρ_lw) - r_lw_soil;

    return nothing
end;

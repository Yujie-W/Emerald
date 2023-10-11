# This file contains function to compute the canopy structure related parameters before running the canopy optical models

#######################################################################################################################################################################################################
#
# Changes to this file
# General
#     2023-Oct-10: add function canopy_structure! (run only once per inclination angle distribution)
#     2023-Oct-11: compute longwave reflectance, transmittance, and emissivity
#
#######################################################################################################################################################################################################
"""

    canopy_structure!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Update canopy structure related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function canopy_structure!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; DIM_LAYER, Θ_INCL) = config;
    (; CANOPY, LEAVES, SOIL_BULK) = spac;

    # compute the weighed average of the leaf inclination angle distribution
    CANOPY.structure.auxil.bf = 0;
    for i in eachindex(Θ_INCL)
        CANOPY.structure.auxil.bf += CANOPY.structure.state.p_incl[i] * cosd(Θ_INCL[i]) ^ 2;;
    end;
    CANOPY.structure.auxil.ddb = (1 + CANOPY.structure.auxil.bf) / 2;
    CANOPY.structure.auxil.ddf = (1 - CANOPY.structure.auxil.bf) / 2;

    # compute longwave effective emissivity, reflectance, and transmittance per layer without correction (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    for i in 1:DIM_LAYER
        ilai = CANOPY.structure.state.δlai[i] * CANOPY.structure.auxil.ci;
        σ_lw_b = CANOPY.structure.auxil.ddb * LEAVES[i].bio.auxil.ρ_lw + CANOPY.structure.auxil.ddf * LEAVES[i].bio.auxil.τ_lw;
        σ_lw_f = CANOPY.structure.auxil.ddf * LEAVES[i].bio.auxil.ρ_lw + CANOPY.structure.auxil.ddb * LEAVES[i].bio.auxil.τ_lw;
        CANOPY.structure.auxil.τ_lw_layer[i] = exp(-1 * (1 - σ_lw_f) * ilai);
        CANOPY.structure.auxil.ρ_lw_layer[i] = 1 - exp(-1 * σ_lw_b * ilai);
        CANOPY.structure.auxil.ϵ_lw_layer[i] = 1 - CANOPY.structure.auxil.τ_lw_layer[i] - CANOPY.structure.auxil.ρ_lw_layer[i];
    end;

    # update the effective longwave reflectance and transmittance
    CANOPY.structure.auxil.ρ_lw[end] = SOIL_BULK.auxil.ρ_lw;
    for i in DIM_LAYER:-1:1
        denom = 1 - CANOPY.structure.auxil.ρ_lw_layer[i] * CANOPY.structure.auxil.ρ_lw[i+1];
        CANOPY.structure.auxil.τ_lw[i] = CANOPY.structure.auxil.τ_lw_layer[i] / denom;                                                                          # it, rescale
        CANOPY.structure.auxil.ρ_lw[i] = CANOPY.structure.auxil.ρ_lw_layer[i] + CANOPY.structure.auxil.τ_lw_layer[i] ^ 2 * CANOPY.structure.auxil.ρ_lw[i+1];    # ir + it-jr-it
    end;

    return nothing
end;

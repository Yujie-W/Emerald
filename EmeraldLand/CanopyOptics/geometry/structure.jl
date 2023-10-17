# This file contains function to compute the canopy structure related parameters before running the canopy optical models

#######################################################################################################################################################################################################
#
# Changes to this file
# General
#     2023-Oct-10: add function canopy_structure! (run only once per inclination angle distribution)
#     2023-Oct-11: compute longwave reflectance, transmittance, and emissivity
#     2023-Oct-14: if LAI <= 0, do nothing
#
#######################################################################################################################################################################################################
"""

    canopy_structure!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update canopy structure related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function canopy_structure!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    if spac.CANOPY.structure.state.lai <= 0
        return nothing
    end;

    # run the canopy structure function only when LAI > 0
    (; DIM_LAYER, Θ_INCL) = config;
    (; CANOPY, LEAVES, SOIL_BULK) = spac;

    # compute the weighed average of the leaf inclination angle distribution
    CANOPY.structure.auxil.bf = 0;
    for i in eachindex(Θ_INCL)
        CANOPY.structure.auxil.bf += CANOPY.structure.state.p_incl[i] * cosd(Θ_INCL[i]) ^ 2;;
    end;
    CANOPY.structure.auxil.ddb = (1 + CANOPY.structure.auxil.bf) / 2;
    CANOPY.structure.auxil.ddf = (1 - CANOPY.structure.auxil.bf) / 2;

    # update the clumping index
    # note here the clumping index should not be a solar zenith angle dependent variable for diffuse light !!!
    # TODO: redesign the clumping index model
    # CANOPY.structure.auxil.ci = CANOPY.structure.state.Ω_A + CANOPY.structure.state.Ω_B * (1 - cosd(CANOPY.sun_geometry.state.sza));
    CANOPY.structure.auxil.ci = CANOPY.structure.state.Ω_A;

    # compute the scattering coefficients for the solar radiation per leaf area
    for i in 1:DIM_LAYER
        j = DIM_LAYER + 1 - i;
        CANOPY.structure.auxil.ddb_leaf[:,i] .= CANOPY.structure.auxil.ddb * LEAVES[j].bio.auxil.ρ_leaf .+ CANOPY.structure.auxil.ddf * LEAVES[j].bio.auxil.τ_leaf;
        CANOPY.structure.auxil.ddf_leaf[:,i] .= CANOPY.structure.auxil.ddf * LEAVES[j].bio.auxil.ρ_leaf .+ CANOPY.structure.auxil.ddb * LEAVES[j].bio.auxil.τ_leaf;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    CANOPY.structure.auxil.τ_dd_layer .= exp.(-1 .* (1 .- CANOPY.structure.auxil.ddf_leaf) .* CANOPY.structure.state.δlai' .* CANOPY.structure.auxil.ci);
    CANOPY.structure.auxil.ρ_dd_layer .= 1 .- exp.(-1 .* CANOPY.structure.auxil.ddb_leaf .* CANOPY.structure.state.δlai' .* CANOPY.structure.auxil.ci);

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    CANOPY.structure.auxil.ρ_dd[:,end] .= SOIL_BULK.auxil.ρ_sw;
    for i in DIM_LAYER:-1:1
        ρ_dd_layer = view(CANOPY.structure.auxil.ρ_dd_layer,:,i  );
        ρ_dd_i     = view(CANOPY.structure.auxil.ρ_dd      ,:,i  );
        ρ_dd_j     = view(CANOPY.structure.auxil.ρ_dd      ,:,i+1);
        τ_dd_layer = view(CANOPY.structure.auxil.τ_dd_layer,:,i  );
        τ_dd_i     = view(CANOPY.structure.auxil.τ_dd      ,:,i  );

        τ_dd_i .= τ_dd_layer ./ (1 .- ρ_dd_layer .* ρ_dd_j);            # ddit; rescale
        ρ_dd_i .= ρ_dd_layer .+ τ_dd_layer .* ρ_dd_j .* τ_dd_i;         # ddir + ddit-ddjr-ddit
    end;

    # compute longwave effective emissivity, reflectance, and transmittance per layer without correction (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    for i in 1:DIM_LAYER
        j = DIM_LAYER + 1 - i;
        ilai = CANOPY.structure.state.δlai[i] * CANOPY.structure.auxil.ci;
        σ_lw_b = CANOPY.structure.auxil.ddb * LEAVES[j].bio.auxil.ρ_lw + CANOPY.structure.auxil.ddf * LEAVES[j].bio.auxil.τ_lw;
        σ_lw_f = CANOPY.structure.auxil.ddf * LEAVES[j].bio.auxil.ρ_lw + CANOPY.structure.auxil.ddb * LEAVES[j].bio.auxil.τ_lw;
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

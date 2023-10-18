# This file contains function to compute the canopy structure related parameters before running the canopy optical models

#######################################################################################################################################################################################################
#
# Changes to this file
# General
#     2023-Oct-10: add function canopy_structure! (run only once per inclination angle distribution)
#     2023-Oct-11: compute longwave reflectance, transmittance, and emissivity
#     2023-Oct-14: if LAI <= 0, do nothing
#     2023-Oct-18: account for SAI in the canopy structure calculation
#
#######################################################################################################################################################################################################
"""

    canopy_structure!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update canopy structure related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function canopy_structure!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;

    if can_str.state.lai <= 0 && can_str.state.sai <= 0
        return nothing
    end;

    # run the canopy structure function only when LAI+SAI > 0
    (; DIM_LAYER, SPECTRA, Θ_INCL) = config;
    sbulk = spac.soil_bulk;

    # compute the weighed average of the leaf inclination angle distribution
    can_str.auxil.bf = 0;
    for i in eachindex(Θ_INCL)
        can_str.auxil.bf += can_str.state.p_incl[i] * cosd(Θ_INCL[i]) ^ 2;;
    end;
    can_str.auxil.ddb = (1 + can_str.auxil.bf) / 2;
    can_str.auxil.ddf = (1 - can_str.auxil.bf) / 2;

    # update the clumping index
    # note here the clumping index should not be a solar zenith angle dependent variable for diffuse light !!!
    # TODO: redesign the clumping index model
    # can_str.auxil.ci = can_str.state.Ω_A + can_str.state.Ω_B * (1 - cosd(spac.canopy.sun_geometry.state.sza));
    can_str.auxil.ci = can_str.state.Ω_A;

    # compute the scattering coefficients for the solar radiation per leaf area
    for i in 1:DIM_LAYER
        leaf = spac.plant.leaves[DIM_LAYER + 1 - i];
        can_str.auxil.ddb_leaf[:,i] .= can_str.auxil.ddb * leaf.bio.auxil.ρ_leaf .+ can_str.auxil.ddf * leaf.bio.auxil.τ_leaf;
        can_str.auxil.ddf_leaf[:,i] .= can_str.auxil.ddf * leaf.bio.auxil.ρ_leaf .+ can_str.auxil.ddb * leaf.bio.auxil.τ_leaf;
        can_str.auxil.ddb_stem[:,i] .= can_str.auxil.ddb * SPECTRA.ρ_STEM;
        can_str.auxil.ddf_stem[:,i] .= can_str.auxil.ddf * SPECTRA.ρ_STEM;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    # can_str.auxil.τ_dd_layer .= exp.(-1 .* (1 .- can_str.auxil.ddf_leaf) .* can_str.state.δlai' .* can_str.auxil.ci);
    # can_str.auxil.ρ_dd_layer .= 1 .- exp.(-1 .* can_str.auxil.ddb_leaf .* can_str.state.δlai' .* can_str.auxil.ci);
    for i in 1:DIM_LAYER
        k_τ_x = (can_str.state.δlai[i] .* (1 .- can_str.auxil.ddf_leaf[:,i]) .+ can_str.state.δsai[i] .* (1 .- can_str.auxil.ddf_stem[:,i])) .* can_str.auxil.ci;
        k_ρ_x = (can_str.state.δlai[i] .* can_str.auxil.ddb_leaf[:,i] .+ can_str.state.δsai[i] .* can_str.auxil.ddb_stem[:,i]) .* can_str.auxil.ci;
        can_str.auxil.τ_dd_layer[:,i] .= exp.(-1 .* k_τ_x);
        can_str.auxil.ρ_dd_layer[:,i] .= 1 .- exp.(-1 .* k_ρ_x);
    end;

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    can_str.auxil.ρ_dd[:,end] .= sbulk.auxil.ρ_sw;
    for i in DIM_LAYER:-1:1
        ρ_dd_layer = view(can_str.auxil.ρ_dd_layer,:,i  );
        ρ_dd_i     = view(can_str.auxil.ρ_dd      ,:,i  );
        ρ_dd_j     = view(can_str.auxil.ρ_dd      ,:,i+1);
        τ_dd_layer = view(can_str.auxil.τ_dd_layer,:,i  );
        τ_dd_i     = view(can_str.auxil.τ_dd      ,:,i  );

        τ_dd_i .= τ_dd_layer ./ (1 .- ρ_dd_layer .* ρ_dd_j);            # ddit; rescale
        ρ_dd_i .= ρ_dd_layer .+ τ_dd_layer .* ρ_dd_j .* τ_dd_i;         # ddir + ddit-ddjr-ddit
    end;

    # compute longwave effective emissivity, reflectance, and transmittance per layer without correction (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    for i in 1:DIM_LAYER
        leaf = spac.plant.leaves[DIM_LAYER + 1 - i];
        ipai = (can_str.state.δlai[i] + can_str.state.δsai[i]) * can_str.auxil.ci;
        σ_lw_b = can_str.auxil.ddb * leaf.bio.auxil.ρ_lw + can_str.auxil.ddf * leaf.bio.auxil.τ_lw;
        σ_lw_f = can_str.auxil.ddf * leaf.bio.auxil.ρ_lw + can_str.auxil.ddb * leaf.bio.auxil.τ_lw;
        can_str.auxil.τ_lw_layer[i] = exp(-1 * (1 - σ_lw_f) * ipai);
        can_str.auxil.ρ_lw_layer[i] = 1 - exp(-1 * σ_lw_b * ipai);
        can_str.auxil.ϵ_lw_layer[i] = 1 - can_str.auxil.τ_lw_layer[i] - can_str.auxil.ρ_lw_layer[i];
    end;

    # update the effective longwave reflectance and transmittance
    can_str.auxil.ρ_lw[end] = sbulk.auxil.ρ_lw;
    for i in DIM_LAYER:-1:1
        denom = 1 - can_str.auxil.ρ_lw_layer[i] * can_str.auxil.ρ_lw[i+1];
        can_str.auxil.τ_lw[i] = can_str.auxil.τ_lw_layer[i] / denom;                                                              # it, rescale
        can_str.auxil.ρ_lw[i] = can_str.auxil.ρ_lw_layer[i] + can_str.auxil.τ_lw_layer[i] ^ 2 * can_str.auxil.ρ_lw[i+1];    # ir + it-jr-it
    end;

    return nothing
end;

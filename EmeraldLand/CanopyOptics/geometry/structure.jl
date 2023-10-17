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
    can_struct = spac.canopy.structure;

    if can_struct.state.lai <= 0
        return nothing
    end;

    # run the canopy structure function only when LAI > 0
    (; DIM_LAYER, Θ_INCL) = config;
    soil_bulk = spac.soil_bulk;

    # compute the weighed average of the leaf inclination angle distribution
    can_struct.auxil.bf = 0;
    for i in eachindex(Θ_INCL)
        can_struct.auxil.bf += can_struct.state.p_incl[i] * cosd(Θ_INCL[i]) ^ 2;;
    end;
    can_struct.auxil.ddb = (1 + can_struct.auxil.bf) / 2;
    can_struct.auxil.ddf = (1 - can_struct.auxil.bf) / 2;

    # update the clumping index
    # note here the clumping index should not be a solar zenith angle dependent variable for diffuse light !!!
    # TODO: redesign the clumping index model
    # can_struct.auxil.ci = can_struct.state.Ω_A + can_struct.state.Ω_B * (1 - cosd(spac.canopy.sun_geometry.state.sza));
    can_struct.auxil.ci = can_struct.state.Ω_A;

    # compute the scattering coefficients for the solar radiation per leaf area
    for i in 1:DIM_LAYER
        leaf = spac.plant.leaves[DIM_LAYER + 1 - i];
        can_struct.auxil.ddb_leaf[:,i] .= can_struct.auxil.ddb * leaf.bio.auxil.ρ_leaf .+ can_struct.auxil.ddf * leaf.bio.auxil.τ_leaf;
        can_struct.auxil.ddf_leaf[:,i] .= can_struct.auxil.ddf * leaf.bio.auxil.ρ_leaf .+ can_struct.auxil.ddb * leaf.bio.auxil.τ_leaf;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    can_struct.auxil.τ_dd_layer .= exp.(-1 .* (1 .- can_struct.auxil.ddf_leaf) .* can_struct.state.δlai' .* can_struct.auxil.ci);
    can_struct.auxil.ρ_dd_layer .= 1 .- exp.(-1 .* can_struct.auxil.ddb_leaf .* can_struct.state.δlai' .* can_struct.auxil.ci);

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    can_struct.auxil.ρ_dd[:,end] .= soil_bulk.auxil.ρ_sw;
    for i in DIM_LAYER:-1:1
        ρ_dd_layer = view(can_struct.auxil.ρ_dd_layer,:,i  );
        ρ_dd_i     = view(can_struct.auxil.ρ_dd      ,:,i  );
        ρ_dd_j     = view(can_struct.auxil.ρ_dd      ,:,i+1);
        τ_dd_layer = view(can_struct.auxil.τ_dd_layer,:,i  );
        τ_dd_i     = view(can_struct.auxil.τ_dd      ,:,i  );

        τ_dd_i .= τ_dd_layer ./ (1 .- ρ_dd_layer .* ρ_dd_j);            # ddit; rescale
        ρ_dd_i .= ρ_dd_layer .+ τ_dd_layer .* ρ_dd_j .* τ_dd_i;         # ddir + ddit-ddjr-ddit
    end;

    # compute longwave effective emissivity, reflectance, and transmittance per layer without correction (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    for i in 1:DIM_LAYER
        leaf = spac.plant.leaves[DIM_LAYER + 1 - i];
        ilai = can_struct.state.δlai[i] * can_struct.auxil.ci;
        σ_lw_b = can_struct.auxil.ddb * leaf.bio.auxil.ρ_lw + can_struct.auxil.ddf * leaf.bio.auxil.τ_lw;
        σ_lw_f = can_struct.auxil.ddf * leaf.bio.auxil.ρ_lw + can_struct.auxil.ddb * leaf.bio.auxil.τ_lw;
        can_struct.auxil.τ_lw_layer[i] = exp(-1 * (1 - σ_lw_f) * ilai);
        can_struct.auxil.ρ_lw_layer[i] = 1 - exp(-1 * σ_lw_b * ilai);
        can_struct.auxil.ϵ_lw_layer[i] = 1 - can_struct.auxil.τ_lw_layer[i] - can_struct.auxil.ρ_lw_layer[i];
    end;

    # update the effective longwave reflectance and transmittance
    can_struct.auxil.ρ_lw[end] = soil_bulk.auxil.ρ_lw;
    for i in DIM_LAYER:-1:1
        denom = 1 - can_struct.auxil.ρ_lw_layer[i] * can_struct.auxil.ρ_lw[i+1];
        can_struct.auxil.τ_lw[i] = can_struct.auxil.τ_lw_layer[i] / denom;                                                              # it, rescale
        can_struct.auxil.ρ_lw[i] = can_struct.auxil.ρ_lw_layer[i] + can_struct.auxil.τ_lw_layer[i] ^ 2 * can_struct.auxil.ρ_lw[i+1];    # ir + it-jr-it
    end;

    return nothing
end;

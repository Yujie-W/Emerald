# This file contains function to compute the sun geometry related parameters (apart from sensor geometry)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sun_geometry! (run per solar zenith angle)
#     2023-Oct-11: compute canopy layer scattering, reflectance, and transmittance
#     2023-Oct-14: do nothing if sza > 89 or LAI <= 0
#     2023-Oct-18: account for SAI in the sun geometry calculation
#
#######################################################################################################################################################################################################
"""

    sun_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update sun geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sun_geometry!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;
    sun_geo = spac.canopy.sun_geometry;

    if sun_geo.state.sza > 89 || (can_str.trait.lai <= 0 && can_str.trait.sai <= 0)
        return nothing
    end;

    # if sza <= 89 and LAI+SAI > 0, run the sun geometry function
    (; SPECTRA) = config;
    leaves = spac.plant.leaves;
    sbulk = spac.soil_bulk;
    n_layer = length(leaves);

    # compute the scattering coefficients for the solar radiation per leaf area
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sun_geo.auxil.sdb_leaf[:,irt] .= sun_geo.s_aux.sdb * leaf.bio.auxil.ρ_leaf .+ sun_geo.s_aux.sdf * leaf.bio.auxil.τ_leaf;
        sun_geo.auxil.sdf_leaf[:,irt] .= sun_geo.s_aux.sdf * leaf.bio.auxil.ρ_leaf .+ sun_geo.s_aux.sdb * leaf.bio.auxil.τ_leaf;
        sun_geo.auxil.sdb_stem[:,irt] .= sun_geo.s_aux.sdb * SPECTRA.ρ_STEM;
        sun_geo.auxil.sdf_stem[:,irt] .= sun_geo.s_aux.sdf * SPECTRA.ρ_STEM;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    # sun_geo.auxil.τ_ss_layer .= exp.(-1 .* sun_geo.s_aux.ks .* can_str.trait.δlai .* can_str.trait.ci);
    # sun_geo.auxil.τ_sd_layer .= 1 .- exp.(-1 .* sun_geo.auxil.sdf_leaf .* can_str.trait.δlai' .* can_str.trait.ci);
    # sun_geo.auxil.ρ_sd_layer .= 1 .- exp.(-1 .* sun_geo.auxil.sdb_leaf .* can_str.trait.δlai' .* can_str.trait.ci);
    for i in 1:n_layer
        kt_ss_x = sun_geo.s_aux.ks .* (can_str.trait.δlai[i] + can_str.trait.δsai[i]) .* can_str.trait.ci;
        kt_sd_x = (sun_geo.auxil.sdf_leaf[:,i] .* can_str.trait.δlai[i] .+ sun_geo.auxil.sdf_stem[:,i] .* can_str.trait.δsai[i]) .* can_str.trait.ci;
        kr_sd_x = (sun_geo.auxil.sdb_leaf[:,i] .* can_str.trait.δlai[i] .+ sun_geo.auxil.sdb_stem[:,i] .* can_str.trait.δsai[i]) .* can_str.trait.ci;
        sun_geo.auxil.τ_ss_layer[i] = exp(-kt_ss_x);
        sun_geo.auxil.τ_sd_layer[:,i] .= 1 .- exp.(-1 .* kt_sd_x);
        sun_geo.auxil.ρ_sd_layer[:,i] .= 1 .- exp.(-1 .* kr_sd_x);
    end;

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    sun_geo.auxil.ρ_sd[:,end] .= sbulk.auxil.ρ_sw;
    for i in n_layer:-1:1
        ρ_dd_layer = view(can_str.auxil.ρ_dd_layer,:,i  );
        ρ_dd_j     = view(can_str.auxil.ρ_dd      ,:,i+1);
        ρ_sd_layer = view(sun_geo.auxil.ρ_sd_layer,:,i  );
        ρ_sd_i     = view(sun_geo.auxil.ρ_sd      ,:,i  );
        ρ_sd_j     = view(sun_geo.auxil.ρ_sd      ,:,i+1);
        τ_dd_layer = view(can_str.auxil.τ_dd_layer,:,i  );
        τ_sd_layer = view(sun_geo.auxil.τ_sd_layer,:,i  );
        τ_sd_i     = view(sun_geo.auxil.τ_sd      ,:,i  );
        τ_ss_layer = view(sun_geo.auxil.τ_ss_layer,  i  );

        τ_sd_i .= (τ_sd_layer .+ τ_ss_layer .* ρ_sd_j .* ρ_dd_layer) ./ (1 .- ρ_dd_layer .* ρ_dd_j);    # sdit + ssit-sdjr-ddit; rescale
        ρ_sd_i .= ρ_sd_layer .+ τ_ss_layer .* ρ_sd_j .* τ_dd_layer .+ τ_sd_i .* ρ_dd_j .* τ_dd_layer;   # sdir + ssit-sdjr-ddit + sdit-ddjr-ddit
    end;

    return nothing
end;

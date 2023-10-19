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

    if sun_geo.state.sza > 89 || (can_str.state.lai <= 0 && can_str.state.sai <= 0)
        return nothing
    end;

    # if sza <= 89 and LAI+SAI > 0, run the sun geometry function
    (; DIM_LAYER, SPECTRA, Θ_AZI, Θ_INCL) = config;
    leaves = spac.plant.leaves;
    sbulk = spac.soil_bulk;

    # extinction coefficients for the solar radiation
    for i in eachindex(Θ_INCL)
        Cs = cosd(Θ_INCL[i]) * cosd(sun_geo.state.sza);
        Ss = sind(Θ_INCL[i]) * sind(sun_geo.state.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        sun_geo.auxil.Cs_incl[i] = Cs;
        sun_geo.auxil.Ss_incl[i] = Ss;
        sun_geo.auxil.βs_incl[i] = βs;
        sun_geo.auxil.ks_incl[i] = 2 / FT(π) / cosd(sun_geo.state.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;
    sun_geo.auxil.ks = can_str.auxil.p_incl' * sun_geo.auxil.ks_incl;

    # compute the scattering weights for diffuse/direct -> diffuse for backward and forward scattering
    sun_geo.auxil.sdb = (sun_geo.auxil.ks + can_str.auxil.bf) / 2;
    sun_geo.auxil.sdf = (sun_geo.auxil.ks - can_str.auxil.bf) / 2;

    # compute the fs and fs_abs matrices
    for i in eachindex(Θ_AZI)
        view(sun_geo.auxil.fs,:,i) .= sun_geo.auxil.Cs_incl .+ sun_geo.auxil.Ss_incl .* cosd(Θ_AZI[i]);
    end;
    sun_geo.auxil.fs ./= cosd(sun_geo.state.sza);
    sun_geo.auxil.fs_abs .= abs.(sun_geo.auxil.fs);
    mul!(sun_geo.auxil.fs_abs_mean, sun_geo.auxil.fs_abs', can_str.auxil.p_incl);
    for i in eachindex(Θ_INCL)
        view(sun_geo.auxil.fs_cos²_incl,i,:) .= view(sun_geo.auxil.fs,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;

    # compute the sunlit leaf fraction
    # sun_geo.auxil.ps = exp.(sun_geo.auxil.ks .* can_str.auxil.ci * can_str.state.lai .* can_str.auxil.x_bnds);
    kscipai = sun_geo.auxil.ks * can_str.auxil.ci * (can_str.state.lai + can_str.state.sai);
    for i in 1:DIM_LAYER
        ksciipai = sun_geo.auxil.ks * can_str.auxil.ci * (can_str.state.δlai[i] + can_str.state.δsai[i]);
        sun_geo.auxil.p_sunlit[i] = 1 / ksciipai * (exp(kscipai * can_str.auxil.x_bnds[i]) - exp(kscipai * can_str.auxil.x_bnds[i+1]));
    end;

    # compute the scattering coefficients for the solar radiation per leaf area
    for i in 1:DIM_LAYER
        leaf = leaves[DIM_LAYER + 1 - i];
        sun_geo.auxil.sdb_leaf[:,i] .= sun_geo.auxil.sdb * leaf.bio.auxil.ρ_leaf .+ sun_geo.auxil.sdf * leaf.bio.auxil.τ_leaf;
        sun_geo.auxil.sdf_leaf[:,i] .= sun_geo.auxil.sdf * leaf.bio.auxil.ρ_leaf .+ sun_geo.auxil.sdb * leaf.bio.auxil.τ_leaf;
        sun_geo.auxil.sdb_stem[:,i] .= sun_geo.auxil.sdb * SPECTRA.ρ_STEM;
        sun_geo.auxil.sdf_stem[:,i] .= sun_geo.auxil.sdf * SPECTRA.ρ_STEM;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    # sun_geo.auxil.τ_ss_layer .= exp.(-1 .* sun_geo.auxil.ks .* can_str.state.δlai .* can_str.auxil.ci);
    # sun_geo.auxil.τ_sd_layer .= 1 .- exp.(-1 .* sun_geo.auxil.sdf_leaf .* can_str.state.δlai' .* can_str.auxil.ci);
    # sun_geo.auxil.ρ_sd_layer .= 1 .- exp.(-1 .* sun_geo.auxil.sdb_leaf .* can_str.state.δlai' .* can_str.auxil.ci);
    for i in 1:DIM_LAYER
        kt_ss_x = sun_geo.auxil.ks .* (can_str.state.δlai[i] + can_str.state.δsai[i]) .* can_str.auxil.ci;
        kt_sd_x = (sun_geo.auxil.sdf_leaf[:,i] .* can_str.state.δlai[i] .+ sun_geo.auxil.sdf_stem[:,i] .* can_str.state.δsai[i]) .* can_str.auxil.ci;
        kr_sd_x = (sun_geo.auxil.sdb_leaf[:,i] .* can_str.state.δlai[i] .+ sun_geo.auxil.sdb_stem[:,i] .* can_str.state.δsai[i]) .* can_str.auxil.ci;
        sun_geo.auxil.τ_ss_layer[i] = exp(-kt_ss_x);
        sun_geo.auxil.τ_sd_layer[:,i] .= 1 .- exp.(-1 .* kt_sd_x);
        sun_geo.auxil.ρ_sd_layer[:,i] .= 1 .- exp.(-1 .* kr_sd_x);
    end;

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    sun_geo.auxil.ρ_sd[:,end] .= sbulk.auxil.ρ_sw;
    for i in DIM_LAYER:-1:1
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

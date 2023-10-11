# This file contains function to compute the sun geometry related parameters (apart from sensor geometry)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sun_geometry! (run per solar zenith angle)
#     2023-Oct-10: compute clumping index as part of the sun geometry function
#     2023-Oct-11: compute canopy layer scattering, reflectance, and transmittance
#
#######################################################################################################################################################################################################
"""

    sun_geometry!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Update sun geometry related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function sun_geometry!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; DIM_LAYER, Θ_AZI, Θ_INCL) = config;
    (; CANOPY, LEAVES, SOIL_BULK) = spac;

    # extinction coefficients for the solar radiation
    for i in eachindex(Θ_INCL)
        Cs = cosd(Θ_INCL[i]) * cosd(CANOPY.sun_geometry.state.sza);
        Ss = sind(Θ_INCL[i]) * sind(CANOPY.sun_geometry.state.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        CANOPY.sun_geometry.auxil.Cs_incl[i] = Cs;
        CANOPY.sun_geometry.auxil.Ss_incl[i] = Ss;
        CANOPY.sun_geometry.auxil.βs_incl[i] = βs;
        CANOPY.sun_geometry.auxil.ks_incl[i] = 2 / FT(π) / cosd(CANOPY.sun_geometry.state.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;
    CANOPY.sun_geometry.auxil.ks = CANOPY.structure.state.p_incl' * CANOPY.sun_geometry.auxil.ks_incl;

    # compute the scattering weights for diffuse/direct -> diffuse for backward and forward scattering
    CANOPY.sun_geometry.auxil.sdb = (CANOPY.sun_geometry.auxil.ks + CANOPY.structure.auxil.bf) / 2;
    CANOPY.sun_geometry.auxil.sdf = (CANOPY.sun_geometry.auxil.ks - CANOPY.structure.auxil.bf) / 2;

    # compute the fs and fs_abs matrices
    for i in eachindex(Θ_AZI)
        view(CANOPY.sun_geometry.auxil.fs,:,i) .= CANOPY.sun_geometry.auxil.Cs_incl .+ CANOPY.sun_geometry.auxil.Ss_incl .* cosd(Θ_AZI[i]);
    end;
    CANOPY.sun_geometry.auxil.fs ./= cosd(CANOPY.sun_geometry.state.sza);
    CANOPY.sun_geometry.auxil.fs_abs .= abs.(CANOPY.sun_geometry.auxil.fs);
    for i in eachindex(Θ_INCL)
        view(CANOPY.sun_geometry.auxil.fs_cos²_incl,i,:) .= view(CANOPY.sun_geometry.auxil.fs,i,:) * cosd(Θ_INCL[i]) ^ 2;
    end;

    # update the clumping index
    CANOPY.structure.auxil.ci = CANOPY.structure.state.Ω_A + CANOPY.structure.state.Ω_B * (1 - cosd(CANOPY.sun_geometry.state.sza));

    # compute the sunlit leaf fraction
    CANOPY.sun_geometry.auxil.ps = exp.(CANOPY.sun_geometry.auxil.ks .* CANOPY.structure.auxil.ci * CANOPY.structure.state.lai .* CANOPY.structure.auxil.x_bnds);
    for i in 1:DIM_LAYER
        ilai = CANOPY.structure.auxil.ci * CANOPY.structure.state.δlai[i];
        kscilai = CANOPY.sun_geometry.auxil.ks * CANOPY.structure.auxil.ci * CANOPY.structure.state.lai;
        CANOPY.sun_geometry.auxil.p_sunlit[i] = 1 / (CANOPY.sun_geometry.auxil.ks * ilai) * (exp(kscilai * CANOPY.structure.auxil.x_bnds[i]) - exp(kscilai * CANOPY.structure.auxil.x_bnds[i+1]));
    end;

    # compute the scattering coefficients for the solar radiation per leaf area
    for i in eachindex(LEAVES)
        CANOPY.sun_geometry.auxil.ddb_leaf[:,i] .= CANOPY.structure.auxil.ddb * LEAVES[i].bio.auxil.ρ_leaf .+ CANOPY.structure.auxil.ddf * LEAVES[i].bio.auxil.τ_leaf;
        CANOPY.sun_geometry.auxil.ddf_leaf[:,i] .= CANOPY.structure.auxil.ddf * LEAVES[i].bio.auxil.ρ_leaf .+ CANOPY.structure.auxil.ddb * LEAVES[i].bio.auxil.τ_leaf;
        CANOPY.sun_geometry.auxil.sdb_leaf[:,i] .= CANOPY.sun_geometry.auxil.sdb * LEAVES[i].bio.auxil.ρ_leaf .+ CANOPY.sun_geometry.auxil.sdf * LEAVES[i].bio.auxil.τ_leaf;
        CANOPY.sun_geometry.auxil.sdf_leaf[:,i] .= CANOPY.sun_geometry.auxil.sdf * LEAVES[i].bio.auxil.ρ_leaf .+ CANOPY.sun_geometry.auxil.sdb * LEAVES[i].bio.auxil.τ_leaf;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    CANOPY.sun_geometry.auxil.τ_ss_layer .= exp.(-1 .* CANOPY.sun_geometry.auxil.ks .* CANOPY.structure.state.δlai .* CANOPY.structure.auxil.ci);
    CANOPY.sun_geometry.auxil.τ_dd_layer .= exp.(-1 .* (1 .- CANOPY.sun_geometry.auxil.ddf_leaf) .* CANOPY.structure.state.δlai' .* CANOPY.structure.auxil.ci);
    CANOPY.sun_geometry.auxil.τ_sd_layer .= 1 .- exp.(-1 .* CANOPY.sun_geometry.auxil.sdf_leaf .* CANOPY.structure.state.δlai' .* CANOPY.structure.auxil.ci);
    CANOPY.sun_geometry.auxil.ρ_dd_layer .= 1 .- exp.(-1 .* CANOPY.sun_geometry.auxil.ddb_leaf .* CANOPY.structure.state.δlai' .* CANOPY.structure.auxil.ci);
    CANOPY.sun_geometry.auxil.ρ_sd_layer .= 1 .- exp.(-1 .* CANOPY.sun_geometry.auxil.sdb_leaf .* CANOPY.structure.state.δlai' .* CANOPY.structure.auxil.ci);

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    CANOPY.sun_geometry.auxil.ρ_dd[:,end] .= SOIL_BULK.auxil.ρ_sw;
    CANOPY.sun_geometry.auxil.ρ_sd[:,end] .= SOIL_BULK.auxil.ρ_sw;
    for i in DIM_LAYER:-1:1
        ρ_dd_layer = view(CANOPY.sun_geometry.auxil.ρ_dd_layer,:,i  );
        ρ_dd_i     = view(CANOPY.sun_geometry.auxil.ρ_dd      ,:,i  );
        ρ_dd_j     = view(CANOPY.sun_geometry.auxil.ρ_dd      ,:,i+1);
        ρ_sd_layer = view(CANOPY.sun_geometry.auxil.ρ_sd_layer,:,i  );
        ρ_sd_i     = view(CANOPY.sun_geometry.auxil.ρ_sd      ,:,i  );
        ρ_sd_j     = view(CANOPY.sun_geometry.auxil.ρ_sd      ,:,i+1);
        τ_dd_layer = view(CANOPY.sun_geometry.auxil.τ_dd_layer,:,i  );
        τ_dd_i     = view(CANOPY.sun_geometry.auxil.τ_dd      ,:,i  );
        τ_sd_layer = view(CANOPY.sun_geometry.auxil.τ_sd_layer,:,i  );
        τ_sd_i     = view(CANOPY.sun_geometry.auxil.τ_sd      ,:,i  );
        τ_ss_layer = view(CANOPY.sun_geometry.auxil.τ_ss_layer,  i  );

        τ_sd_i .= (τ_sd_layer .+ τ_ss_layer .* ρ_sd_j .* ρ_dd_layer) ./ (1 .- ρ_dd_layer .* ρ_dd_j);    # sdit + ssit-sdjr-ddit; rescale
        τ_dd_i .= τ_dd_layer ./ (1 .- ρ_dd_layer .* ρ_dd_j);                                            # ddit; rescale
        ρ_sd_i .= ρ_sd_layer .+ τ_ss_layer .* ρ_sd_j .* τ_dd_layer .+ τ_sd_i .* ρ_dd_j .* τ_dd_layer;   # sdir + ssit-sdjr-ddit + sdit-ddjr-ddit
        ρ_dd_i .= ρ_dd_layer .+ τ_dd_layer .* ρ_dd_j .* τ_dd_i;                                         # ddir + ddit-ddjr-ddit
    end;

    return nothing
end;

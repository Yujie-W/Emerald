# This file contains function to compute the canopy structure related parameters before running the canopy optical models

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function canopy_structure_aux! to update the trait-dependent auxiliary variables for canopy structure (to call in t_aux!)
#     2024-Mar-01: compute the extinction coefficients for the diffuse radiation (isotropic)
#     2024-Mar-01: compute the fraction of ddb and ddf (relative fraction of purely backward and forward scattering)
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: compute CI weighted extinction coefficients for the diffuse radiation (so no need to multiply by CI in the equations when using kd_leaf and kd_stem)
#     2024-Sep-09: compute the CI for the diffuse radiation based on the angle dependent CI
#     2024-Oct-16: compute the diffuse light extinction coefficient per 1 degree (so no need to weigh it here)
#
#######################################################################################################################################################################################################
"""

    canopy_structure_aux!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}

Update the trait-dependent auxiliary variables for canopy structure, given
- `config` SPAC configuration
- `can` MultiLayerCanopy

"""
function canopy_structure_aux! end;

canopy_structure_aux!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT} = canopy_structure_aux!(config, can.structure.trait, can.structure.t_aux);

canopy_structure_aux!(config::SPACConfiguration{FT}, trait::CanopyStructureTrait{FT}, t_aux::CanopyStructureTDAuxil{FT}) where {FT} = (
    (; Θ_INCL, Θ_INCL_BNDS) = config;

    # update the clumping index
    # note here the clumping index should not be a solar zenith angle dependent variable for diffuse light
    # thus, ci is computed as the integral of zenith angle dependent ci of the diffuse light
    #     ci = ∫_0^π/2 ci(θ) * sind(θ) dθ / ∫_0^π/2 sind(θ) dθ
    #     ∫_0^π/2 ci(θ) * sind(θ) dθ = ci_0 * ∫_0^π/2 sind(θ) dθ - ci_0 * ci_1 * ∫_0^π/2 cos(θ) sind(θ) dθ = ci_0 - ci_0 * ci_1 / 2
    #     ∫_0^π/2 sind(θ) dθ = 1
    #     ci = ci_0 - ci_0 * ci_1 / 2
    t_aux.ci_diffuse = trait.ci.ci_0 - trait.ci.ci_0 * trait.ci.ci_1 / 2;

    # compute the probability of leaf inclination angles based on lidf
    for i in eachindex(t_aux.p_incl_leaf)
        t_aux.p_incl_leaf[i] = lidf_cdf(trait.lidf, Θ_INCL_BNDS[i,2]) - lidf_cdf(trait.lidf, Θ_INCL_BNDS[i,1]);
        t_aux.p_incl_stem[i] = lidf_cdf(trait.sidf, Θ_INCL_BNDS[i,2]) - lidf_cdf(trait.sidf, Θ_INCL_BNDS[i,1]);
    end;

    # compute the extinction coefficients for the diffuse radiation (directions in isotropic radiation)
    # these kd_leaf and kd_stem already account for the impact of clumping index
    for i_dif in 1:90
        θ_dif = i_dif - FT(0.5);
        kd_l = 0;
        kd_s = 0;
        for i in eachindex(Θ_INCL)
            kd_l += extinction_coefficient(θ_dif, Θ_INCL[i], trait.ci) * t_aux.p_incl_leaf[i];
            kd_s += extinction_coefficient(θ_dif, Θ_INCL[i], trait.ci) * t_aux.p_incl_stem[i];
        end;
        t_aux.kd_leaf[i_dif] = kd_l;
        t_aux.kd_stem[i_dif] = kd_s;
    end;

    # compute the weighed average of the leaf inclination angle distribution
    t_aux.ddb_leaf = 0;
    t_aux.ddf_leaf = 0;
    t_aux.ddb_stem = 0;
    t_aux.ddf_stem = 0;
    for i in eachindex(Θ_INCL)
        f_ada = f_adaxial(Θ_INCL[i]);
        f_aba = 1 - f_ada;
        f_inc = Θ_INCL[i] / 180;
        t_aux.ddb_leaf += (f_ada * (1 - f_inc) + f_aba * f_inc) * t_aux.p_incl_leaf[i];
        t_aux.ddf_leaf += (f_ada * f_inc + f_aba * (1 - f_inc)) * t_aux.p_incl_leaf[i];
        t_aux.ddb_stem += (f_ada * (1 - f_inc) + f_aba * f_inc) * t_aux.p_incl_stem[i];
        t_aux.ddf_stem += (f_ada * f_inc + f_aba * (1 - f_inc)) * t_aux.p_incl_stem[i];
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this file
# General
#     2023-Oct-10: add function canopy_structure! (run only once per inclination angle distribution)
#     2023-Oct-11: compute longwave reflectance, transmittance, and emissivity
#     2023-Oct-14: if LAI <= 0, do nothing
#     2023-Oct-18: account for SAI in the canopy structure calculation
#     2024-Mar-01: compute the layer shortwave and longwave scattering coefficients based on the new theory
#     2024-Sep-04: separate leaf and stem optical properties
#     2024-Sep-07: compute CI weighted extinction coefficients for the diffuse radiation (so no need to multiply by CI in the equations when using kd_leaf and kd_stem)
#     2024-Oct-16: weigh the extinction coefficient for diffuse radiation when computing the transmittance
#     2024-Oct-16: add option to compute effective leaf spectra based on CI
#     2024-Nov-08: when using EFFECTIVE_LEAF_SPECTRA make sure LAI > 0
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

    if can_str.trait.lai <= 0 && can_str.trait.sai <= 0
        return nothing
    end;

    # run the canopy structure function only when LAI+SAI > 0
    (; EFFECTIVE_LEAF_SPECTRA, SPECTRA) = config;
    leaves = spac.plant.leaves;
    sbulk = spac.soil_bulk;
    n_layer = length(leaves);

    # compute the effective leaf reflectance and transmittance spectra using the PROSPECT model scheme
    # Note: this step can be bypassed if leaf reflectance and transmittance are not changing
    mask_effective = EFFECTIVE_LEAF_SPECTRA && can_str.trait.lai > 0;
    if mask_effective
        ρ_2 = spac.cache.cache_wl_1;
        τ_2 = spac.cache.cache_wl_2;
        n_eff = 1 / can_str.t_aux.ci_diffuse;
        for irt in 1:n_layer
            ilf = n_layer + 1 - irt;
            leaf = leaves[ilf];
            ρ_1 = leaf.bio.auxil.ρ_leaf;
            τ_1 = leaf.bio.auxil.τ_leaf;
            ρ_2 .= layer_2_ρ.(ρ_1, τ_1, n_eff - 1);
            τ_2 .= layer_2_τ.(ρ_1, τ_1, n_eff - 1);
            can_str.auxil.ρ_leaf_eff[:,irt] .= leaf_ρ.(ρ_1, τ_1, ρ_1, τ_1, ρ_2);
            can_str.auxil.τ_leaf_eff[:,irt] .= leaf_τ.(τ_1, ρ_1, ρ_2, τ_2);
        end;
    end;

    # compute the scattering coefficients for the solar radiation per leaf area
    # TODO: add cosine of the diffuse angle and CI impacts on the rho and tau of an effective leaf
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = spac.plant.leaves[ilf];
        ρ_leaf = mask_effective ? view(can_str.auxil.ρ_leaf_eff,:,irt) : leaf.bio.auxil.ρ_leaf;
        τ_leaf = mask_effective ? view(can_str.auxil.τ_leaf_eff,:,irt) : leaf.bio.auxil.τ_leaf;
        can_str.auxil.ddb_leaf[:,irt] .= can_str.t_aux.ddb_leaf .* ρ_leaf .+ can_str.t_aux.ddf_leaf .* τ_leaf;
        can_str.auxil.ddf_leaf[:,irt] .= can_str.t_aux.ddf_leaf .* ρ_leaf .+ can_str.t_aux.ddb_leaf .* τ_leaf;
        can_str.auxil.ddb_stem[:,irt] .= can_str.t_aux.ddb_stem .* SPECTRA.ρ_STEM;
        can_str.auxil.ddf_stem[:,irt] .= can_str.t_aux.ddf_stem .* SPECTRA.ρ_STEM;
    end;

    # compute the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    #     can_str.auxil.τ_dd_layer .= exp.(-1 .* (1 .- can_str.auxil.ddf_leaf) .* can_str.trait.δlai' .* can_str.trait.ci);
    #     can_str.auxil.ρ_dd_layer .= 1 .- exp.(-1 .* can_str.auxil.ddb_leaf .* can_str.trait.δlai' .* can_str.trait.ci);
    # Later, we included support to SAI as well.
    # However, as of 2024-Feb-29, we found an issue with the equations above when LAI and SAI are big enough in a single layer that sum of reflectance and transmittance is greater than 1.
    # Therefore, we revised the equations using calculus and the equations for the whole layer are as follows (for the transmitted and reflected light):
    #     can_str.auxil.ρ_dd_layer = ∫_0^iCIPAI (ddb_leaf * δLAI + ddb_stem * δSAI) / δPAI * kd * exp(-kd * x) * dx = (ddb_leaf * δLAI + ddb_stem * δSAI) / δPAI * (1 - exp(-kd * iCIPAI))
    #     can_str.auxil.τ_dd_layer = ∫_0^iCIPAI (ddf_leaf * δLAI + ddf_stem * δSAI) / δPAI * kd * exp(-kd * x) * dx = (ddf_leaf * δLAI + ddf_stem * δSAI) / δPAI * (1 - exp(-kd * iCIPAI))
    # Then, the total transmitted radiation need to plus the radiation that has not passed though any leaf, namely τ_dd_solar.
    # As of 2024-Sep-07, we account CI's angular dependency along with kd_leaf and kd_stem, and thus CI is removed here from the equations below.
    # As of 2024-Oct-16, we realized that the extinction coefficient for diffuse radiation (weighted with sind(θ_za)) may overestimate the extinction coefficient for the diffuse radiation.
    # Therefore, we compute the extinction coefficient per 1 degree and weigh the final transmittance (that do not reach any leaf surface).
    k_τ_x = spac.cache.cache_wl_1;
    k_ρ_x = spac.cache.cache_wl_2;
    for irt in 1:n_layer
        δlai = can_str.trait.δlai[irt];
        δsai = can_str.trait.δsai[irt];
        δpai = δlai + δsai;
        # loop through the isotropic directions
        τ_dd_weighed = 0;
        sum_sind = 0;
        for i_dif in 1:90
            θ_dif = i_dif - FT(0.5);
            kt_dd_i_dif = can_str.t_aux.kd_leaf[i_dif] * δlai + can_str.t_aux.kd_stem[i_dif] * δsai;
            τ_dd_i_dif = exp(-kt_dd_i_dif);
            τ_dd_weighed += τ_dd_i_dif * sind(θ_dif);
            sum_sind += sind(θ_dif);
        end;
        τ_dd_weighed /= sum_sind;
        can_str.auxil.τ_dd_isotropic[irt] = τ_dd_weighed;
        k_τ_x .= (view(can_str.auxil.ddf_leaf,:,irt) .* δlai .+ view(can_str.auxil.ddf_stem,:,irt) .* δsai) ./ δpai;
        k_ρ_x .= (view(can_str.auxil.ddb_leaf,:,irt) .* δlai .+ view(can_str.auxil.ddb_stem,:,irt) .* δsai) ./ δpai;
        can_str.auxil.τ_dd_layer[:,irt] .= (1 - τ_dd_weighed) .* k_τ_x .+ τ_dd_weighed;
        can_str.auxil.ρ_dd_layer[:,irt] .= (1 - τ_dd_weighed) .* k_ρ_x;
    end;

    # compute the effective tranmittance and reflectance per layer from lowest to highest layer (including the denominator correction)
    can_str.auxil.ρ_dd[:,end] .= sbulk.auxil.ρ_sw;
    for i in n_layer:-1:1
        ρ_dd_layer = view(can_str.auxil.ρ_dd_layer,:,i  );
        ρ_dd_i     = view(can_str.auxil.ρ_dd      ,:,i  );
        ρ_dd_j     = view(can_str.auxil.ρ_dd      ,:,i+1);
        τ_dd_layer = view(can_str.auxil.τ_dd_layer,:,i  );
        τ_dd_i     = view(can_str.auxil.τ_dd      ,:,i  );

        τ_dd_i .= τ_dd_layer ./ (1 .- ρ_dd_layer .* ρ_dd_j);        # ddit; rescale
        ρ_dd_i .= ρ_dd_layer .+ τ_dd_layer .* ρ_dd_j .* τ_dd_i;     # ddir + ddit-ddjr-ddit
    end;

    # compute longwave effective emissivity, reflectance, and transmittance per layer without correction (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    # Similar to the bug fix related to diffuse radiation, we fix the equations for the longwave radiation as well.
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = spac.plant.leaves[ilf];
        δlai = can_str.trait.δlai[irt];
        δsai = can_str.trait.δsai[irt];
        δpai = δlai + δsai;
        σ_leaf_b = can_str.t_aux.ddb_leaf * leaf.bio.trait.ρ_lw + can_str.t_aux.ddf_leaf * leaf.bio.trait.τ_lw;
        σ_leaf_f = can_str.t_aux.ddf_leaf * leaf.bio.trait.ρ_lw + can_str.t_aux.ddb_leaf * leaf.bio.trait.τ_lw;
        σ_stem_b = can_str.t_aux.ddb_stem * leaf.bio.trait.ρ_lw;
        σ_stem_f = can_str.t_aux.ddf_stem * leaf.bio.trait.ρ_lw;
        k_ρ_x = (σ_leaf_b * δlai .+ σ_stem_b * δsai) ./ δpai;
        k_τ_x = (σ_leaf_f * δlai .+ σ_stem_f * δsai) ./ δpai;
        τ_dd_lw = exp(-can_str.auxil.τ_dd_isotropic[irt]);
        can_str.auxil.τ_lw_layer[irt] = (1 - τ_dd_lw) .* k_τ_x .+ τ_dd_lw;
        can_str.auxil.ρ_lw_layer[irt] = (1 - τ_dd_lw) .* k_ρ_x;
        can_str.auxil.ϵ_lw_layer[irt] = 1 - can_str.auxil.τ_lw_layer[irt] - can_str.auxil.ρ_lw_layer[irt];
    end;

    # update the effective longwave reflectance and transmittance
    can_str.auxil.ρ_lw[end] = sbulk.trait.ρ_lw;
    for i in n_layer:-1:1
        denom = 1 - can_str.auxil.ρ_lw_layer[i] * can_str.auxil.ρ_lw[i+1];
        can_str.auxil.τ_lw[i] = can_str.auxil.τ_lw_layer[i] / denom;                                                        # it, rescale
        can_str.auxil.ρ_lw[i] = can_str.auxil.ρ_lw_layer[i] + can_str.auxil.τ_lw_layer[i] ^ 2 * can_str.auxil.ρ_lw[i+1];    # ir + it-jr-it
    end;

    return nothing
end;

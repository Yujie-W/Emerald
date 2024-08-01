#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function to photosystem_coefficients!
#     2022-Jan-14: add function that operates PSM, PRC, and FLM directly so as to be more modular (reduce memory allocations)
#     2022-Feb-07: use ppar in fluorescence model
#     2022-Feb-07: add support for Johnson and Berry (2021) model
#     2022-Feb-07: use a_gross and j_pot rather than a series of j_p680 and j_p700
#     2022-Feb-10: scale fluorescence quantum yield based on f_psii and reabsorption factor
#     2022-Feb-10: q1 needs to be multiply by η
#     2022-Mar-04: add support to sustained NPQ
#     2022-Mar-04: use the weighted yield for photosynthesis
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Jun-15: set ϕ_f and ϕ_p to 0 when ppar is 0
#     2023-Oct-24: save PSI and PSII ϕ_f in the C3Cyto model
#     2023-Oct-28: add method for QLFluoscenceModel
#     2023-Oct-30: compute q_l using exp(-flm.K_B * ppar) (omitting K_A)
#     2024-Aug-01: generalize the function for GeneralC3Trait and GeneralC4Trait
#     2024-Aug-01: compute q_l using flm.K_A * exp(-flm.K_B * ppar) for QLFluoscenceModelHan model
# Bug fixes
#     2022-Feb-24: a typo from "rc.ϕ_f  = rc.f_m′ / (1 - rc.ϕ_p);" to "rc.ϕ_f  = rc.f_m′ * (1 - rc.ϕ_p);"
#     2022-Feb-28: ps.e_to_c is recalculated based on analytically resolving leaf.p_CO₂_i from leaf.g_CO₂, this ps.e_to_c used to be calculated as ps.a_j / ps.j (a_j here is not p_CO₂_i based)
#                  note here that in CliMA v0.1, this e_to_c is not updated properly.
#     2024-Feb-26: use k_d from auxil rather than from trait (because k_d is updated in the auxil to account for its TD)
# To do
#     TODO: add more calculations such as NPQ when the model is ready
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(ps::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, ppar::FT; β::FT = FT(1)) where {FT}

Update the rate constants and coefficients in reaction center, given
- `ps` `C3Cyto`, `C3VJP`, or `C4VJP` type photosynthesis system
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_coefficients! end;

# Leaf
photosystem_coefficients!(
            config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(config, ps.trait, ps.state, ps.auxil, ppar; β = β);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(config, pss, psa, pst.FLM, ppar; β = β);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            flm::CytochromeFluoscenceModel{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PSI_RATE_CONSTANTS, PSII_RATE_CONSTANTS) = config;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    ϕ_P1_a = psa.a_g * psa.η / (psa.e2c * ppar * (1 - psa.f_psii));
    ϕ_P2_a = psa.a_g / (psa.e2c * ppar * psa.f_psii);
    q1     = ϕ_P1_a / psa.ϕ_psi_max;
    q2     = 1 - psa.j_psi / (β * psa.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    k_sum_nb = -1 * (PSII_RATE_CONSTANTS.K_U * ϕ_P2_a + PSII_RATE_CONSTANTS.K_P * (q2 - ϕ_P2_a));
    k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * PSII_RATE_CONSTANTS.K_U * PSII_RATE_CONSTANTS.K_P);
    k_sum    = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    k_n      = k_sum - PSII_RATE_CONSTANTS.K_F - PSII_RATE_CONSTANTS.K_U - PSII_RATE_CONSTANTS.K_D;

    # compute PSII and PSI yeilds
    k_sum_1 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + k_n;
    k_sum_2 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + k_n + PSII_RATE_CONSTANTS.K_P;
    k_sum_3 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P;
    k_sum_4 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_X;
    ϕ_U2_a  = (q2 * PSII_RATE_CONSTANTS.K_U / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_U / k_sum_1);
    ϕ_F2_a  = (q2 * PSII_RATE_CONSTANTS.K_F / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_F / k_sum_1) / (1 - ϕ_U2_a);
    ϕ_F1_a  = PSI_RATE_CONSTANTS.K_F / k_sum_3 * q1 + PSI_RATE_CONSTANTS.K_F / k_sum_4 * (1 - q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    psa.ϕ_f1 = ϕ_F1_a;
    psa.ϕ_f2 = ϕ_F2_a;
    psa.ϕ_f  = ϕ_F1_a * (1 - psa.f_psii) + ϕ_F2_a * psa.f_psii;
    psa.ϕ_p  = ϕ_P1_a * (1 - psa.f_psii) + ϕ_P2_a * psa.f_psii;

    return nothing
);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            flm::KNFluoscenceModel{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = config;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate rate constants
    x  = max(0, 1 - psa.ϕ_p / psa.ϕ_psii_max);
    xᵅ = x ^ flm.K_A;
    psa.k_n = flm.K_0 * (1 + flm.K_B) * xᵅ / (flm.K_B + xᵅ);
    psa.k_p = max(0, psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus) / (1 - psa.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P + psa.k_d);
    psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            flm::QLFluoscenceModel{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = config;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = exp(-flm.K_B * ppar);
    psa.k_p = PS_RATE_CONSTANTS.K_P * q_l;
    psa.k_n = (psa.k_p - psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P + psa.k_d);
    psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            flm::QLFluoscenceModelHan{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = config;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = flm.K_A * exp(-flm.K_B * ppar);
    psa.k_p = PS_RATE_CONSTANTS.K_P * q_l;
    psa.k_n = (psa.k_p - psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P + psa.k_d);
    psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

# For CanopyLayer
photosystem_coefficients!(
            config::SPACConfiguration{FT},
            cache::SPACCache{FT},
            ps::CanopyLayerPhotosystem{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(config, cache, ps.trait, ps.state, ps.auxil, ppar; β = β);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            cache::SPACCache{FT},
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(config, cache, pss, psa, pst.FLM, ppar; β = β);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            cache::SPACCache{FT},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            flm::CytochromeFluoscenceModel{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    (; PSI_RATE_CONSTANTS, PSII_RATE_CONSTANTS) = config;

    # unpack the vars from the cache
    ϕ_P1_a   = cache.cache_incl_azi_2_1;
    ϕ_P2_a   = cache.cache_incl_azi_2_2;
    q1       = cache.cache_incl_azi_2_3;
    q2       = cache.cache_incl_azi_2_4;
    k_sum_nb = cache.cache_incl_azi_2_5;
    k_sum_nc = cache.cache_incl_azi_2_6;
    k_sum    = cache.cache_incl_azi_2_7;
    k_n      = cache.cache_incl_azi_2_8;
    k_sum_1  = cache.cache_incl_azi_2_9;
    k_sum_2  = cache.cache_incl_azi_3_1;
    ϕ_U2_a   = cache.cache_incl_azi_3_2;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    @. ϕ_P1_a = psa.a_g * psa.η / (psa.e2c * ppar * (1 - psa.f_psii));
    @. ϕ_P2_a = psa.a_g / (psa.e2c * ppar * psa.f_psii);
    @. q1     = ϕ_P1_a / psa.ϕ_psi_max;
    @. q2     = 1 - psa.j_psi / (β * psa.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    @. k_sum_nb = -1 * (PSII_RATE_CONSTANTS.K_U * ϕ_P2_a + PSII_RATE_CONSTANTS.K_P * (q2 - ϕ_P2_a));
    @. k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * PSII_RATE_CONSTANTS.K_U * PSII_RATE_CONSTANTS.K_P);
    @. k_sum    = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    @. k_n      = k_sum - PSII_RATE_CONSTANTS.K_F - PSII_RATE_CONSTANTS.K_U - PSII_RATE_CONSTANTS.K_D;

    # compute PSII and PSI yeilds
    @. k_sum_1 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + k_n;
    @. k_sum_2 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + k_n + PSII_RATE_CONSTANTS.K_P;
    k_sum_3 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P;
    k_sum_4 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_X;
    @. ϕ_U2_a = (q2 * PSII_RATE_CONSTANTS.K_U / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_U / k_sum_1);
    @. psa.ϕ_f2 = (q2 * PSII_RATE_CONSTANTS.K_F / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_F / k_sum_1) / (1 - ϕ_U2_a);
    @. psa.ϕ_f1 = PSI_RATE_CONSTANTS.K_F / k_sum_3 * q1 + PSI_RATE_CONSTANTS.K_F / k_sum_4 * (1 - q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    @. psa.ϕ_f = psa.ϕ_f1 * (1 - psa.f_psii) + psa.ϕ_f2 * psa.f_psii;
    @. psa.ϕ_p = ϕ_P1_a * (1 - psa.f_psii) + ϕ_P2_a * psa.f_psii;

    return nothing
);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            cache::SPACCache{FT},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            flm::KNFluoscenceModel{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = config;

    # calculate photochemical yield
    psa.ϕ_p .= psa.a_g ./ (psa.e2c .* psa.f_psii .* ppar);

    # unpack the cache
    x = cache.cache_incl_azi_2_1;
    xᵅ = cache.cache_incl_azi_2_2;

    # calculate rate constants
    @. x  = max(0, 1 - psa.ϕ_p / psa.ϕ_psii_max);
    @. xᵅ = x ^ flm.K_A;
    @. psa.k_n = flm.K_0 .* (1 + flm.K_B) .* xᵅ ./ (flm.K_B .+ xᵅ);
    @. psa.k_p = max.(0, psa.ϕ_p .* (PS_RATE_CONSTANTS.K_F .+ psa.k_d .+ psa.k_n .+ pss.k_npq_sus) ./ (1 .- psa.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    psa.f_o = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + PS_RATE_CONSTANTS.K_P);
    psa.f_m = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    @. psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus + PS_RATE_CONSTANTS.K_P);
    @. psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            cache::SPACCache{FT},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            flm::QLFluoscenceModel{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = config;

    # calculate photochemical yield
    @. psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = cache.cache_incl_azi_2_1;
    @. q_l = exp(-flm.K_B * ppar);
    @. psa.k_p = PS_RATE_CONSTANTS.K_P * q_l;
    @. psa.k_n = (psa.k_p - psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + PS_RATE_CONSTANTS.K_P);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    @. psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus + PS_RATE_CONSTANTS.K_P);
    @. psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            config::SPACConfiguration{FT},
            cache::SPACCache{FT},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            flm::QLFluoscenceModelHan{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = config;

    # calculate photochemical yield
    @. psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = cache.cache_incl_azi_2_1;
    @. q_l = flm.K_A * exp(-flm.K_B * ppar);
    @. psa.k_p = PS_RATE_CONSTANTS.K_P * q_l;
    @. psa.k_n = (psa.k_p - psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + PS_RATE_CONSTANTS.K_P);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    @. psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus + PS_RATE_CONSTANTS.K_P);
    @. psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

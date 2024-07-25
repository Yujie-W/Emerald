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
#     2023-Sep-09: compute ϕ_d and ϕ_n in the VJPReactionCenter
#     2023-Oct-24: save PSI and PSII ϕ_f in the C3Cyto model
#     2023-Oct-28: add method for QLFluoscenceModel
#     2023-Oct-30: compute q_l using exp(-flm.K_B * ppar) (omitting K_A)
#     2024-Jul-23: add support to C3CytoMinEtaTrait
# Bug fixes
#     2022-Feb-24: a typo from "rc.ϕ_f  = rc.f_m′ / (1 - rc.ϕ_p);" to "rc.ϕ_f  = rc.f_m′ * (1 - rc.ϕ_p);"
#     2022-Feb-28: psm.e_to_c is recalculated based on analytically resolving leaf.p_CO₂_i from leaf.g_CO₂, this psm.e_to_c used to be calculated as psm.a_j / psm.j (a_j here is not p_CO₂_i based)
#                  note here that in CliMA v0.1, this e_to_c is not updated properly.
#     2024-Feb-26: use k_d from auxil rather than from trait (because k_d is updated in the auxil to account for its TD)
# To do
#     TODO: add more calculations such as NPQ when the model is ready
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(psm::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, ppar::FT; β::FT = FT(1)) where {FT}

Update the rate constants and coefficients in reaction center, given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` type photosynthesis system
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_coefficients! end;

photosystem_coefficients!(psm::LeafPhotosystem{FT}, ppar::FT; β::FT = FT(1)) where {FT} = photosystem_coefficients!(psm.trait, psm.state, psm.auxil, ppar; β = β);

photosystem_coefficients!(pst::Union{C3CytoMinEtaTrait{FT}, C3CytoTrait{FT}, C3JBTrait{FT}}, pss::C3State{FT}, psa::LeafPhotosystemAuxil{FT}, ppar::FT; β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    ϕ_P1_a = psa.a_g * psa.η / (psa.e2c * ppar * (1 - psa.f_psii));
    ϕ_P2_a = psa.a_g / (psa.e2c * ppar * psa.f_psii);
    q1     = ϕ_P1_a / psa.ϕ_psi_max;
    q2     = 1 - psa.j_psi / (β * psa.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    k_sum_nb = -1 * (pst.K_U * ϕ_P2_a + pst.K_PSII * (q2 - ϕ_P2_a));
    k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * pst.K_U * pst.K_PSII);
    k_sum    = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    k_n      = k_sum - pst.K_F - pst.K_U - pst.K_D;

    # compute PSII and PSI yeilds
    k_sum_1 = pst.K_D + pst.K_F + pst.K_U + k_n;
    k_sum_2 = pst.K_D + pst.K_F + pst.K_U + k_n + pst.K_PSII;
    k_sum_3 = pst.K_D + pst.K_F + pst.K_PSI;
    k_sum_4 = pst.K_D + pst.K_F + pst.K_X;
    ϕ_U2_a  =  q2 * pst.K_U / k_sum_2 + (1 - q2) * pst.K_U / k_sum_1;
    ϕ_F2_a  = (q2 * pst.K_F / k_sum_2 + (1 - q2) * pst.K_F / k_sum_1) / (1 - ϕ_U2_a);
    ϕ_F1_a  = pst.K_F / k_sum_3 * q1 + pst.K_F / k_sum_4 * (1 - q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    psa.ϕ_f1 = ϕ_F1_a;
    psa.ϕ_f2 = ϕ_F2_a;
    psa.ϕ_f  = ϕ_F1_a * (1 - psa.f_psii) + ϕ_F2_a * psa.f_psii;
    psa.ϕ_p  = ϕ_P1_a * (1 - psa.f_psii) + ϕ_P2_a * psa.f_psii;

    return nothing
);

photosystem_coefficients!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} =
    photosystem_coefficients!(pst, pss, pst.FLM, psa, ppar; β = β);

photosystem_coefficients!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            flm::KNFluoscenceModel{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate rate constants
    x  = max(0, 1 - psa.ϕ_p / psa.ϕ_psii_max);
    xᵅ = x ^ flm.K_A;
    psa.k_n = flm.K_0 * (1 + flm.K_B) * xᵅ / (flm.K_B + xᵅ);
    psa.k_p = max(0, psa.ϕ_p * (pst.K_F + psa.k_d + psa.k_n + pss.k_npq_sus) / (1 - psa.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    psa.f_o  = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d);
    psa.f_o′ = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.f_m  = pst.K_F / (pst.K_F + psa.k_d);
    psa.f_m′ = pst.K_F / (pst.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_d  = psa.k_d / pst.K_F * psa.ϕ_f;
    psa.ϕ_n  = (psa.k_n + pss.k_npq_sus) / pst.K_F * psa.ϕ_f;
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_n + pss.k_npq_sus) / (pst.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            flm::QLFluoscenceModel{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT;
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = exp(-flm.K_B * ppar);
    psa.k_p = pst.K_PSII * q_l;
    psa.k_n = (psa.k_p - psa.ϕ_p * (pst.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d);
    psa.f_o′ = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.f_m  = pst.K_F / (pst.K_F + psa.k_d);
    psa.f_m′ = pst.K_F / (pst.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_d  = psa.k_d / pst.K_F * psa.ϕ_f;
    psa.ϕ_n  = (psa.k_n + pss.k_npq_sus) / pst.K_F * psa.ϕ_f;
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_n + pss.k_npq_sus) / (pst.K_F + psa.k_d);

    return nothing
);

# For CanopyLayer
photosystem_coefficients!(psm::CanopyLayerPhotosystem{FT}, ppar::Vector{FT}; β::FT = FT(1)) where {FT} = photosystem_coefficients!(psm.trait, psm.state, psm.auxil, ppar; β = β);

photosystem_coefficients!(pst::Union{C3CytoMinEtaTrait{FT}, C3CytoTrait{FT}, C3JBTrait{FT}}, pss::C3State{FT}, psa::CanopyLayerPhotosystemAuxil{FT}, ppar::Vector{FT}; β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    ϕ_P1_a = psa.a_g * psa.η / (psa.e2c * ppar * (1 - psa.f_psii));
    ϕ_P2_a = psa.a_g / (psa.e2c * ppar * psa.f_psii);
    q1     = ϕ_P1_a / psa.ϕ_psi_max;
    q2     = 1 - psa.j_psi / (β * psa.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    k_sum_nb = -1 * (pst.K_U * ϕ_P2_a + pst.K_PSII * (q2 - ϕ_P2_a));
    k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * pst.K_U * pst.K_PSII);
    k_sum    = upper_quadratic.(k_sum_na, k_sum_nb, k_sum_nc);
    k_n      = k_sum - pst.K_F - pst.K_U - pst.K_D;

    # compute PSII and PSI yeilds
    k_sum_1 = pst.K_D + pst.K_F + pst.K_U + k_n;
    k_sum_2 = pst.K_D + pst.K_F + pst.K_U + k_n + pst.K_PSII;
    k_sum_3 = pst.K_D + pst.K_F + pst.K_PSI;
    k_sum_4 = pst.K_D + pst.K_F + pst.K_X;
    ϕ_U2_a  =  q2 * pst.K_U / k_sum_2 + (1 - q2) * pst.K_U / k_sum_1;
    ϕ_F2_a  = (q2 * pst.K_F / k_sum_2 + (1 - q2) * pst.K_F / k_sum_1) / (1 - ϕ_U2_a);
    ϕ_F1_a  = pst.K_F / k_sum_3 * q1 + pst.K_F / k_sum_4 * (1 - q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    psa.ϕ_f1 = ϕ_F1_a;
    psa.ϕ_f2 = ϕ_F2_a;
    psa.ϕ_f  = ϕ_F1_a * (1 - psa.f_psii) + ϕ_F2_a * psa.f_psii;
    psa.ϕ_p  = ϕ_P1_a * (1 - psa.f_psii) + ϕ_P2_a * psa.f_psii;

    return nothing
);

photosystem_coefficients!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} =
    photosystem_coefficients!(pst, pss, pst.FLM, psa, ppar; β = β);

photosystem_coefficients!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            flm::KNFluoscenceModel{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # calculate photochemical yield
    psa.ϕ_p .= psa.a_g ./ (psa.e2c .* psa.f_psii .* ppar);

    # calculate rate constants
    x  = max.(0, 1 .- psa.ϕ_p ./ psa.ϕ_psii_max);
    xᵅ = x .^ flm.K_A;
    psa.k_n = flm.K_0 .* (1 + flm.K_B) .* xᵅ ./ (flm.K_B .+ xᵅ);
    psa.k_p = max.(0, psa.ϕ_p .* (pst.K_F .+ psa.k_d .+ psa.k_n .+ pss.k_npq_sus) ./ (1 .- psa.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    @. psa.f_o  = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d);
    @. psa.f_o′ = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.f_m  = pst.K_F / (pst.K_F + psa.k_d);
    @. psa.f_m′ = pst.K_F / (pst.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_d  = psa.k_d / pst.K_F * psa.ϕ_f;
    @. psa.ϕ_n  = (psa.k_n + pss.k_npq_sus) / pst.K_F * psa.ϕ_f;
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (pst.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            flm::QLFluoscenceModel{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = exp.(-flm.K_B * ppar);
    psa.k_p = pst.K_PSII * q_l;
    psa.k_n = (psa.k_p - psa.ϕ_p * (pst.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d);
    psa.f_o′ = pst.K_F / (pst.K_F + pst.K_PSII + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.f_m  = pst.K_F / (pst.K_F + psa.k_d);
    psa.f_m′ = pst.K_F / (pst.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_d  = psa.k_d / pst.K_F * psa.ϕ_f;
    psa.ϕ_n  = (psa.k_n + pss.k_npq_sus) / pst.K_F * psa.ϕ_f;
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_n + pss.k_npq_sus) / (pst.K_F + psa.k_d);

    return nothing
);

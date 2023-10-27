#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function to photosystem_coefficients!
#     2022-Jan-14: add function that operates PSM, PRC, and FLM directly so as to be more modular (reduce memory allocations)
#     2022-Feb-07: use ppar in fluorescence model
#     2022-Feb-07: add support for Johnson and Berry (2021) model
#     2022-Feb-07: use a_gross and j_pot rather than a series of j_p680 and j_p700
#     2022-Feb-10: scale fluorescence quantum yield based on F_PSII and reabsorption factor
#     2022-Feb-10: q1 needs to be multiply by η
#     2022-Mar-04: add support to sustained NPQ
#     2022-Mar-04: use the weighted yield for photosynthesis
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Jun-15: set ϕ_f and ϕ_p to 0 when ppar is 0
#     2023-Sep-09: compute ϕ_d and ϕ_n in the VJPReactionCenter
#     2023-Oct-24: save PSI and PSII ϕ_f in the C3Cyto model
# Bug fixes
#     2022-Feb-24: a typo from "rc.ϕ_f  = rc.f_m′ / (1 - rc.ϕ_p);" to "rc.ϕ_f  = rc.f_m′ * (1 - rc.ϕ_p);"
#     2022-Feb-28: psm.e_to_c is recalculated based on analytically resolving leaf.p_CO₂_i from leaf.g_CO₂, this psm.e_to_c used to be calculated as psm.a_j / psm.j (a_j here is not p_CO₂_i based)
#                  note here that in CliMA v0.1, this e_to_c is not updated properly.
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

photosystem_coefficients!(psm::LeafPhotosystem{FT}, ppar::FT; β::FT = FT(1)) where {FT} = photosystem_coefficients!(psm.state, psm.auxil, ppar; β = β);

photosystem_coefficients!(pss::C3CytoState{FT}, psa::PSMAuxil{FT}, ppar::FT; β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    ϕ_P1_a = psa.a_g * psa.η / (psa.e2c * ppar * (1 - pss.F_PSII));
    ϕ_P2_a = psa.a_g / (psa.e2c * ppar * pss.F_PSII);
    q1     = ϕ_P1_a / psa.ϕ_psi_max;
    q2     = 1 - psa.j_psi / (β * psa.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    k_sum_nb = -1 * (pss.K_U * ϕ_P2_a + pss.K_PSII * (q2 - ϕ_P2_a));
    k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * pss.K_U * pss.K_PSII);
    k_sum    = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    k_n      = k_sum - pss.K_F - pss.K_U - pss.K_D;

    # compute PSII and PSI yeilds
    k_sum_1 = pss.K_D + pss.K_F + pss.K_U + k_n;
    k_sum_2 = pss.K_D + pss.K_F + pss.K_U + k_n + pss.K_PSII;
    k_sum_3 = pss.K_D + pss.K_F + pss.K_PSI;
    k_sum_4 = pss.K_D + pss.K_F + pss.K_X;
    ϕ_U2_a  =  q2 * pss.K_U / k_sum_2 + (1 - q2) * pss.K_U / k_sum_1;
    ϕ_F2_a  = (q2 * pss.K_F / k_sum_2 + (1 - q2) * pss.K_F / k_sum_1) / (1 - ϕ_U2_a);
    ϕ_F1_a  = pss.K_F / k_sum_3 * q1 + pss.K_F / k_sum_4 * (1 - q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    psa.ϕ_f1 = ϕ_F1_a;
    psa.ϕ_f2 = ϕ_F2_a;
    psa.ϕ_f  = ϕ_F1_a * (1 - pss.F_PSII) + ϕ_F2_a * pss.F_PSII;
    psa.ϕ_p  = ϕ_P1_a * (1 - pss.F_PSII) + ϕ_P2_a * pss.F_PSII;

    return nothing

    #=
    # some unused or unsaved variables
    _ϕ_N2_a = (_q2 * _k_n  / _k_sum_2 + (1 - _q2) * _k_n / _k_sum_1) / (1 - _ϕ_U2_a);
    _ϕ_D2_a = (_q2 * K_D   / _k_sum_2 + (1 - _q2) * K_D  / _k_sum_1) / (1 - _ϕ_U2_a);
    _ϕ_P1_a = K_PSI / _k_sum_3 * _q1;
    _ϕ_N1_a = K_X   / _k_sum_4 * (1 - _q1);
    _ϕ_D1_a = K_D   / _k_sum_3 * _q1 + K_D / _k_sum_4 * (1 - _q1);

    # PAM measured fluorescence levels (Eqns. 38-42)
    #   N.B., hardcoding of a2(1) for dark-adapted value
    _tmp_1 = α_1 * K_F * ϵ_1;
    _tmp_2 = α_2 * K_F * ϵ_2;
    _Fm_a  = _tmp_1 / _k_sum_4 + _tmp_2 / (K_D + K_F);
    _Fo_a  = _tmp_1 / _k_sum_3 + _tmp_2 / (K_D + K_F + K_PSII);
    _Fmp_a = _tmp_1 / _k_sum_4 + _tmp_2 / (K_D + K_F + _k_n);
    _Fop_a = _tmp_1 / _k_sum_3 + _tmp_2 / (K_D + K_F + K_PSII + _k_n);
    _Fs_a  = α_1 * _ϕ_F1_a * ϵ_1 + α_2 * _ϕ_F2_a * ϵ_2;

    # PAM indices used in plotter_forward_fun.m
    _PAM1_a = 1 - _Fs_a / _Fmp_a; # ϕ_P
    _PAM2_a = _Fs_a * (1 / _Fmp_a - 1/_Fm_a); # ϕ_N
    _PAM3_a = _Fs_a / _Fm_a; # ϕ_D + ϕ_F
    =#
);

photosystem_coefficients!(pss::Union{C3VJPState{FT}, C4VJPState{FT}}, psa::PSMAuxil{FT}, ppar::FT; β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psa.ϕ_f = 0;
        psa.ϕ_p = 0;

        return nothing
    end;

    # calculate photochemical yield
    psa.ϕ_p = psa.a_g / (psa.e2c * pss.F_PSII * ppar);

    # calculate rate constants
    x  = max(0, 1 - psa.ϕ_p / psa.ϕ_psii_max);
    xᵅ = x ^ pss.FLM.K_A;
    psa.k_npq_rev = pss.FLM.K_0 * (1 + pss.FLM.K_B) * xᵅ / (pss.FLM.K_B + xᵅ);
    psa.k_p       = max(0, psa.ϕ_p * (pss.K_F + psa.k_d + psa.k_npq_rev + pss.k_npq_sus) / (1 - psa.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    psa.f_o  = pss.K_F / (pss.K_F + pss.K_PSII + psa.k_d);
    psa.f_o′ = pss.K_F / (pss.K_F + pss.K_PSII + psa.k_d + psa.k_npq_rev + pss.k_npq_sus);
    psa.f_m  = pss.K_F / (pss.K_F + psa.k_d);
    psa.f_m′ = pss.K_F / (pss.K_F + psa.k_d + psa.k_npq_rev + pss.k_npq_sus);
    psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    psa.ϕ_d  = psa.k_d / pss.K_F * psa.ϕ_f;
    psa.ϕ_n  = (psa.k_npq_rev + pss.k_npq_sus) / pss.K_F * psa.ϕ_f;
    psa.ϕ_f1 = psa.ϕ_f;
    psa.ϕ_f2 = psa.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    psa.npq = (psa.k_npq_rev + pss.k_npq_sus) / (pss.K_F + psa.k_d);

    return nothing
);

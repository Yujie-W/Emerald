#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function to photosystem_coefficients!
#     2022-Jan-14: add function that operates PSM, PRC, and FLM directly so as to be more modular (reduce memory allocations)
#     2022-Feb-07: use ppar in fluorescence model
#     2022-Feb-07: add support for Johnson and Berry (2021) model
#     2022-Feb-07: use a_gross and j_pot rather than a series of j_p680 and j_p700
#     2022-Feb-10: scale fluorescence quantum yield based on F_PSI and reabsorption factor
#     2022-Feb-10: q1 needs to be multiply by η
#     2022-Mar-04: add support to sustained NPQ
#     2022-Mar-04: use the weighted yield for photosynthesis
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Jun-15: set ϕ_f and ϕ_p to 0 when ppar is 0
#     2023-Sep-09: compute ϕ_d and ϕ_n in the VJPReactionCenter
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

photosystem_coefficients!(psm::C3Cyto{FT}, ppar::FT; β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psm.auxil.ϕ_f = 0;
        psm.auxil.ϕ_p = 0;

        return nothing
    end;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    ϕ_P1_a = psm.auxil.a_g * psm.auxil.η / (psm.auxil.e2c * ppar * psm.state.F_PSI);
    ϕ_P2_a = psm.auxil.a_g / (psm.auxil.e2c * ppar * (1 - psm.state.F_PSI));
    q1     = ϕ_P1_a / psm.state.Φ_PSI_MAX;
    q2     = 1 - psm.auxil.j_psi / (β * psm.auxil.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    k_sum_nb = -1 * (psm.state.K_U * ϕ_P2_a + psm.state.K_PSII * (q2 - ϕ_P2_a));
    k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * psm.state.K_U * psm.state.K_PSII);
    k_sum    = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    k_n      = k_sum - psm.state.K_F - psm.state.K_U - psm.state.K_D;

    # compute PSII and PSI yeilds
    k_sum_1 = psm.state.K_D + psm.state.K_F + psm.state.K_U + k_n;
    k_sum_2 = psm.state.K_D + psm.state.K_F + psm.state.K_U + k_n + psm.state.K_PSII;
    k_sum_3 = psm.state.K_D + psm.state.K_F + psm.state.K_PSI;
    k_sum_4 = psm.state.K_D + psm.state.K_F + psm.state.K_X;
    ϕ_U2_a  =  q2 * psm.state.K_U / k_sum_2 + (1 - q2) * psm.state.K_U / k_sum_1;
    ϕ_F2_a  = (q2 * psm.state.K_F / k_sum_2 + (1 - q2) * psm.state.K_F / k_sum_1) / (1 - ϕ_U2_a);
    ϕ_F1_a  = psm.state.K_F / k_sum_3 * q1 + psm.state.K_F / k_sum_4 * (1 - q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    psm.auxil.ϕ_f = ϕ_F1_a * psm.auxil.ϵ_1 * psm.state.F_PSI + ϕ_F2_a * psm.auxil.ϵ_2 * (1 - psm.state.F_PSI);
    psm.auxil.ϕ_p = ϕ_P1_a * psm.state.F_PSI + ϕ_P2_a * (1 - psm.state.F_PSI);

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

photosystem_coefficients!(psm::Union{C3VJP{FT}, C4VJP{FT}}, ppar::FT; β::FT = FT(1)) where {FT} = (
    if ppar == 0
        psm.auxil.ϕ_f = 0;
        psm.auxil.ϕ_p = 0;

        return nothing
    end;

    # calculate photochemical yield
    psm.auxil.ϕ_p = psm.auxil.a_g / (psm.auxil.e2c * psm.state.F_PSII * ppar);

    # calculate rate constants
    x  = max(0, 1 - psm.auxil.ϕ_p / psm.auxil.ϕ_psii_max);
    xᵅ = x ^ psm.state.FLM.K_A;
    psm.auxil.k_npq_rev = psm.state.FLM.K_0 * (1 + psm.state.FLM.K_B) * xᵅ / (psm.state.FLM.K_B + xᵅ);
    psm.auxil.k_p       = max(0, psm.auxil.ϕ_p * (psm.state.K_F + psm.auxil.k_d + psm.auxil.k_npq_rev + psm.state.k_npq_sus) / (1 - psm.auxil.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_P_MAX + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_P_MAX + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    psm.auxil.f_o  = psm.state.K_F / (psm.state.K_F + psm.state.K_P_MAX + psm.auxil.k_d);
    psm.auxil.f_o′ = psm.state.K_F / (psm.state.K_F + psm.state.K_P_MAX + psm.auxil.k_d + psm.auxil.k_npq_rev + psm.state.k_npq_sus);
    psm.auxil.f_m  = psm.state.K_F / (psm.state.K_F + psm.auxil.k_d);
    psm.auxil.f_m′ = psm.state.K_F / (psm.state.K_F + psm.auxil.k_d + psm.auxil.k_npq_rev + psm.state.k_npq_sus);
    psm.auxil.ϕ_f  = psm.auxil.f_m′ * (1 - psm.auxil.ϕ_p);
    psm.auxil.ϕ_d  = psm.auxil.k_d / psm.state.K_F * psm.auxil.ϕ_f;
    psm.auxil.ϕ_n  = (psm.auxil.k_npq_rev + psm.state.k_npq_sus) / psm.state.K_F * psm.auxil.ϕ_f;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    psm.auxil.q_e = 1 - (psm.auxil.f_m - psm.auxil.f_o′) / (psm.auxil.f_m′ - psm.auxil.f_o);
    psm.auxil.q_p = 1 - (psm.auxil.ϕ_f - psm.auxil.f_o′) / (psm.auxil.f_m - psm.auxil.f_o′);
    psm.auxil.npq = (psm.auxil.k_npq_rev + psm.state.k_npq_sus) / (psm.state.K_F + psm.auxil.k_d);

    return nothing
);

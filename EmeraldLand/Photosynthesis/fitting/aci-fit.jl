#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to fit A-Ci curve
#     2024-Jul-22: add support to C3CLM, C3FvCB, and C3VJP
#     2024-Jul-23: add support to C3CytoInfApTrait
#     2024-Jul-27: fit Γ_star as well for C3 models
#     2024-Jul-27: add initial guess to option
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#     2024-Aug-06: make initial_guess mandatory
#     2024-Oct-03: customize the initial guess for C3 models
#
#######################################################################################################################################################################################################
"""

    aci_fit(config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            initial_guess::Vector) where {FT}

Fit the A-Ci curve (will be abstractized based on the trait and methods embedded), given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `initial_guess` Initial guess of fitting parameters

"""
function aci_fit end;

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::::Union{Nothing, Vector}) where {FT} = aci_fit(config, ps, ps.trait, air, df, initial_guess);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::::Union{Nothing, Vector}) where {FT} = aci_fit(config, ps, pst, pst.ACM, pst.AJM, pst.APM, air, df, initial_guess);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC3Trait{FT},
        acm::AcMethodC3VcmaxPi,
        ajm::AjMethodC3JmaxPi,
        apm::ApMethodC3Vcmax,
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    if !isnothing(initial_guess) && length(initial_guess) == 4
        mthd = ReduceStepMethodND{FT}(
            x_mins = [1, 1, 1, 0.1],
            x_maxs = [200, 400, 10, 10],
            x_inis = initial_guess,
            Δ_inis = [10, 10, 1, 1],
        );
    else
        # loop through the data once to get the respiration and Γ_star limits
        rd_lim_min = 0.1;
        γ_lim_max = 40;
        for dfr in eachrow(df)
            # estimate the minimum respiration rate limit
            rd_lim_min = max(rd_lim_min, -dfr.A_NET / temperature_correction(ps.trait.TD_R, dfr.T_LEAF));
            # estimate the maximum Γ_star limit
            γ_lim_max = min(γ_lim_max, dfr.P_I / temperature_correction(ps.trait.TD_Γ, dfr.T_LEAF));
        end;
        # loop through the data once again to guess the Vcmax and Jmax
        vcmax_guess = 0;
        jmax_guess = 0;
        ps.trait.r_d25 = rd_lim_min * 1.2;
        ps.trait.TD_Γ.VAL_REF = γ_lim_max * 0.8;
        for dfr in eachrow(df)
            photosystem_temperature_dependence!(config, ps, air, dfr.T_LEAF);
            vcmax_guess = max(vcmax_guess, (dfr.A_NET + ps.auxil.r_d) * (dfr.P_I + ps.auxil.k_m) / (dfr.P_I - ps.auxil.γ_star));
            jmax_guess = max(jmax_guess, (dfr.A_NET + ps.auxil.r_d) * (4*dfr.P_I + 8*ps.auxil.γ_star) / (dfr.P_I - ps.auxil.γ_star) * 1.2);
        end;
        # set the initial guess
        mthd = ReduceStepMethodND{FT}(
            x_mins = [vcmax_guess * 0.5, jmax_guess * 0.5, γ_lim_max * 0.5, rd_lim_min * 1.0],
            x_maxs = [vcmax_guess * 2.0, jmax_guess * 2.0, γ_lim_max * 1.0, rd_lim_min * 2.0],
            x_inis = [vcmax_guess * 1.0, jmax_guess * 1.0, γ_lim_max * 0.8, rd_lim_min * 1.2],
            Δ_inis = [10, 10, 1, 1],
        );
    end;

    stol = SolutionToleranceND{FT}([0.1, 0.1, 0.01, 0.01], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol);
    aci = aci_curve(config, ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC3Trait{FT},
        acm::AcMethodC3VcmaxPi,
        ajm::AjMethodC3VqmaxPi,
        apm::ApMethodC3Vcmax,
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    if !isnothing(initial_guess) && length(initial_guess) == 4
        mthd = ReduceStepMethodND{FT}(
            x_mins = [1, 0.01, 1, 0.1],
            x_maxs = [200, 5, 10, 10],
            x_inis = initial_guess[:],
            Δ_inis = [10, 1, 1, 1],
        );
    else
        # loop through the data once to get the respiration and Γ_star limits
        rd_lim_min = 0.1;
        γ_lim_max = 40;
        for dfr in eachrow(df)
            # estimate the minimum respiration rate limit
            rd_lim_min = max(rd_lim_min, -dfr.A_NET / temperature_correction(ps.trait.TD_R, dfr.T_LEAF));
            # estimate the maximum Γ_star limit
            γ_lim_max = min(γ_lim_max, dfr.P_I / temperature_correction(ps.trait.TD_Γ, dfr.T_LEAF));
        end;
        # loop through the data once again to guess the Vcmax and b6f
        vcmax_guess = 0;
        b6f_guess = 0;
        ps.trait.r_d25 = rd_lim_min * 1.2;
        ps.trait.TD_Γ.VAL_REF = γ_lim_max * 0.8;
        for dfr in eachrow(df)
            photosystem_temperature_dependence!(config, ps, air, dfr.T_LEAF);
            vcmax_guess = max(vcmax_guess, (dfr.A_NET + ps.auxil.r_d) * (dfr.P_I + ps.auxil.k_m) / (dfr.P_I - ps.auxil.γ_star));
            j_1 = (dfr.A_NET + ps.auxil.r_d) * (4*dfr.P_I + 8*ps.auxil.γ_star) / (dfr.P_I - ps.auxil.γ_star) * ps.auxil.η;
            p_1 = dfr.PPAR * 0.5 * ps.auxil.ϕ_psi_max;
            v_q = p_1 * j_1 / (p_1 - j_1);
            b6f_guess = max(b6f_guess, v_q / ps.auxil.k_q);
        end;
        # set the initial guess
        mthd = ReduceStepMethodND{FT}(
            x_mins = [vcmax_guess * 0.5, b6f_guess * 0.5, γ_lim_max * 0.5, rd_lim_min * 1.0],
            x_maxs = [vcmax_guess * 2.0, b6f_guess * 2.0, γ_lim_max * 1.0, rd_lim_min * 2.0],
            x_inis = [vcmax_guess * 1.0, b6f_guess * 1.0, γ_lim_max * 0.8, rd_lim_min * 1.2],
            Δ_inis = [10, 0.1, 1, 1],
        );
    end;
    stol = SolutionToleranceND{FT}([0.1, 0.01, 0.01, 0.01], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol);
    aci = aci_curve(config, ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC4Trait{FT},
        acm::AcMethodC4Vcmax,
        ajm::AjMethodC4JPSII,
        apm::ApMethodC4VcmaxPi,
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::::Union{Nothing, Vector}) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 0.1],
        x_maxs = [200, 10],
        x_inis = initial_guess,
        Δ_inis = [10, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.01], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol);
    aci = aci_curve(config, ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC4Trait{FT},
        acm::AcMethodC4Vcmax,
        ajm::AjMethodC4JPSII,
        apm::ApMethodC4VpmaxPi,
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::::Union{Nothing, Vector}) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 1, 0.1],
        x_maxs = [200, 200, 10],
        x_inis = initial_guess,
        Δ_inis = [10, 10, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.1, 0.01], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol);
    aci = aci_curve(config, ps, air, df);

    return sol, best_rmse, aci
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-20: add functions to fit A-Ci curve with removing outliers
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#     2024-Aug-06: make initial_guess mandatory
#
#######################################################################################################################################################################################################
"""

    aci_fit_exclude_outliter(
                config::SPACConfiguration{FT},
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                df::DataFrame,
                initial_guess::::Union{Nothing, Vector};
                min_count::Int = 9,
                rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve by removing outliers, given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `initial_guess` Initial guess of fitting parameters
- `min_count` Minimum number of data points to fit an A-Ci curve
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit_exclude_outliter(
            config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            initial_guess::::Union{Nothing, Vector};
            min_count::Int = 9,
            rmse_threshold::Number = 2) where {FT}
    # remove outliers using thresholds when necessary
    df[!,"A_NET_BAK"] .= df.A_NET;
    last_rmse = 1000;
    last_sol = nothing;
    last_df = deepcopy(df);
    crnt_df = deepcopy(df);
    while true
        sol, best_rmse, aci = aci_fit(config, ps, air, crnt_df, initial_guess);
        if last_rmse - best_rmse < rmse_threshold
            break
        else
            last_df = deepcopy(crnt_df);
            last_sol = sol;
            last_rmse = best_rmse;
            remove_i = max_diff_index(crnt_df.A_NET, aci);
            if sum(.!isnan.(crnt_df.A_NET)) > min_count
                crnt_df[remove_i, "A_NET"] = NaN;
            end;
        end;
    end;

    # change the df and traits
    df.A_NET .= last_df.A_NET;
    best_rmse = aci_rmse(config, ps, air, df, last_sol);
    aci = aci_curve(config, ps, air, df);

    return last_sol, best_rmse, aci
end;


"""

    max_diff_index(y::Vector, p_pred::Vector)

Find the index of the maximum difference between `y` and `p_pred`.

"""
function max_diff_index(y::Vector, p_pred::Vector)
    abs_diff = abs.(y .- p_pred);
    abs_diff[isnan.(abs_diff)] .= 0;

    return findmax(abs_diff)[2]
end;

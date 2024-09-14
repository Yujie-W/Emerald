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
#
#######################################################################################################################################################################################################
"""

    aci_fit(config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            initial_guess::Vector;
            fit_rd::Bool = false,
            fitted_γ::Union{Nothing, Number} = nothing) where {FT}

Fit the A-Ci curve (will be abstractized based on the trait and methods embedded), given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `initial_guess` Initial guess of fitting parameters
- `fit_rd` Fit Rd or not (last fitting parameter; default: false)
- `fitted_γ` Fitted Γ* value, will be prescribed if given (default: nothing)

"""
function aci_fit end;

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::Vector;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing) where {FT} = aci_fit(config, ps, ps.trait, air, df, initial_guess; fit_rd = fit_rd, fitted_γ = fitted_γ);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::Vector;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing) where {FT} = aci_fit(config, ps, pst, pst.ACM, pst.AJM, pst.APM, air, df, initial_guess; fit_rd = fit_rd, fitted_γ = fitted_γ);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC3Trait{FT},
        acm::AcMethodC3VcmaxPi,
        ajm::AjMethodC3JmaxPi,
        apm::ApMethodC3Vcmax,
        air::AirLayer{FT},
        df::DataFrame,
        initial_guess::Vector;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 1, 1, 0.1] : [1, 1, 1],
        x_maxs = fit_rd ? [200, 400, 10, 10] : [200, 400, 10],
        x_inis = fit_rd ? initial_guess : initial_guess[1:3],
        Δ_inis = fit_rd ? [10, 10, 1, 1] : [10, 10, 1],
    );
    stol = fit_rd ? SolutionToleranceND{FT}([0.1, 0.1, 0.01, 0.01], 50) : SolutionToleranceND{FT}([0.1, 0.1, 0.01], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x; fit_rd = fit_rd, fitted_γ = fitted_γ);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol; fit_rd = fit_rd, fitted_γ = fitted_γ);
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
        initial_guess::Vector;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 0.01, 1, 0.1] : [1, 0.01, 1],
        x_maxs = fit_rd ? [200, 10, 10, 10] : [200, 10, 10],
        x_inis = fit_rd ? initial_guess[:] : initial_guess[1:3],
        Δ_inis = fit_rd ? [10, 1, 1, 1] : [10, 1, 1],
    );
    stol = fit_rd ? SolutionToleranceND{FT}([0.1, 0.01, 0.01, 0.01], 50) : SolutionToleranceND{FT}([0.1, 0.01, 0.01], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x; fit_rd = fit_rd, fitted_γ = fitted_γ);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol; fit_rd = fit_rd, fitted_γ = fitted_γ);
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
        initial_guess::Vector;
        fit_rd::Bool = false,
        fitted_γ::Nothing = nothing) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 0.1] : [1],
        x_maxs = fit_rd ? [200, 10] : [200],
        x_inis = fit_rd ? initial_guess : initial_guess[1:1],
        Δ_inis = fit_rd ? [10, 1] : [10],
    );
    stol = fit_rd ? SolutionToleranceND{FT}([0.1, 0.01], 50) : SolutionToleranceND{FT}([0.1], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x; fit_rd = fit_rd, fitted_γ = fitted_γ);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol; fit_rd = fit_rd, fitted_γ = fitted_γ);
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
        initial_guess::Vector;
        fit_rd::Bool = false,
        fitted_γ::Nothing = nothing) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 1, 0.1] : [1, 1],
        x_maxs = fit_rd ? [200, 200, 10] : [200, 200],
        x_inis = fit_rd ? initial_guess : initial_guess[1:2],
        Δ_inis = fit_rd ? [10, 10, 1] : [10, 10],
    );
    stol = fit_rd ? SolutionToleranceND{FT}([0.1, 0.1, 0.01], 50) : SolutionToleranceND{FT}([0.1, 0.1], 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, x; fit_rd = fit_rd, fitted_γ = fitted_γ);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, sol; fit_rd = fit_rd, fitted_γ = fitted_γ);
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
                initial_guess::Vector;
                fit_rd::Bool = false,
                fitted_γ::Union{Nothing, Number} = nothing,
                min_count::Int = 9,
                rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve by removing outliers, given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `initial_guess` Initial guess of fitting parameters
- `fit_rd` Fit Rd or not (last fitting parameter; default: false)
- `fitted_γ` Fitted Γ* value, will be prescribed if given (default: nothing)
- `min_count` Minimum number of data points to fit an A-Ci curve
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit_exclude_outliter(
            config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            initial_guess::Vector;
            fit_rd::Bool = false,
            fitted_γ::Union{Nothing, Number} = nothing,
            min_count::Int = 9,
            rmse_threshold::Number = 2) where {FT}
    # remove outliers using thresholds when necessary
    df[!,"A_NET_BAK"] .= df.A_NET;
    last_rmse = 1000;
    last_sol = nothing;
    last_df = deepcopy(df);
    crnt_df = deepcopy(df);
    while true
        sol, best_rmse, aci = aci_fit(config, ps, air, crnt_df, initial_guess; fit_rd = fit_rd, fitted_γ = fitted_γ);
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
    best_rmse = aci_rmse(config, ps, air, df, last_sol; fit_rd = fit_rd, fitted_γ = fitted_γ);
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

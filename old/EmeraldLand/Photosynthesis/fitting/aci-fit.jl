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
#     2024-Oct-04: use max(1, p_1 - j_1) to avoid negative b6f guess
#     2025-Jul-30: improved curve fitting algorithm
#
#######################################################################################################################################################################################################
"""

    aci_fit(config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            params::Vector{String},
            initial_guess::Vector) where {FT}

Fit the A-Ci curve (will be abstractized based on the trait and methods embedded), given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `initial_guess` Initial guess of fitting parameters

"""
function aci_fit end;

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = aci_fit(config, ps, ps.trait, air, df, params, initial_guess);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = aci_fit(config, ps, pst, pst.ACM, pst.AJM, pst.APM, air, df, params, initial_guess);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC3Trait{FT},
        acm::AcMethodC3VcmaxPi,
        ajm::AjMethodC3JmaxPi,
        apm::ApMethodC3Vcmax,
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    x_mins = FT[];
    x_maxs = FT[];
    x_inis = FT[];
    Δ_inis = FT[];
    Δ_tols = FT[];
    if !isnothing(initial_guess)
        @assert length(initial_guess) == length(params);
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vcmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Jmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 400);
            iguess = findfirst(params .== "Jmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Γstar25" in params
            push!(x_mins, 1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Γstar25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Rd25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    else
        # loop through the data once to get the respiration and Γ_star limits
        rd_lim_min = 0.1;
        γ_lim_max = 40;
        for dfr in eachrow(df)
            # estimate the minimum respiration rate limit
            rd_lim_min = nanmax([rd_lim_min, -dfr.A_NET / temperature_correction(ps.trait.TD_R, dfr.T_LEAF)]);
            # estimate the maximum Γ_star limit
            γ_lim_max = nanmin([γ_lim_max, dfr.P_I / temperature_correction(ps.trait.TD_Γ, dfr.T_LEAF)]);
        end;
        # loop through the data once again to guess the Vcmax and Jmax
        vcmax_guess = 5;
        jmax_guess = 10;
        ps.trait.r_d25 = rd_lim_min * 1.2;
        ps.trait.TD_Γ.VAL_REF = γ_lim_max * 0.8;
        for dfr in eachrow(df)
            photosystem_temperature_dependence!(config, ps, air, dfr.T_LEAF);
            vcmax = (dfr.A_NET + ps.auxil.r_d) * (dfr.P_I + ps.auxil.k_m) / (dfr.P_I - ps.auxil.γ_star) / temperature_correction(ps.trait.TD_VCMAX, dfr.T_LEAF);
            vcmax_guess = nanmax([vcmax_guess, vcmax]);
            vcmax_guess = nanmin([vcmax_guess, 100]);
            jmax = (dfr.A_NET + ps.auxil.r_d) * (4*dfr.P_I + 8*ps.auxil.γ_star) / (dfr.P_I - ps.auxil.γ_star) * 1.2 / temperature_correction(ps.trait.TD_JMAX, dfr.T_LEAF);
            jmax_guess = nanmax([jmax_guess, jmax]);
            jmax_guess = nanmin([jmax_guess, 200]);
        end;
        # set the initial guess
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, vcmax_guess * 2.0);
            push!(x_inis, vcmax_guess * 1.0);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Jmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, jmax_guess * 2.0);
            push!(x_inis, jmax_guess * 1.0);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Γstar25" in params
            push!(x_mins, 1);
            push!(x_maxs, γ_lim_max * 1.0);
            push!(x_inis, γ_lim_max * 0.8);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Rd25" in params
            push!(x_mins, rd_lim_min * 1.0);
            push!(x_maxs, 5);
            push!(x_inis, rd_lim_min * 1.2);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    end;

    mthd = ReduceStepMethodND{FT}(x_mins = x_mins, x_maxs = x_maxs, x_inis = x_inis, Δ_inis = Δ_inis);
    stol = SolutionToleranceND{FT}(Δ_tols, 50);
    # func(x) = (rme = aci_rmse(config, ps, pst, air, df, x); @info "C3VJP model" x rme ; -rme);
    func(x) = -aci_rmse(config, ps, pst, air, df, params, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, params, sol);
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
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    x_mins = FT[];
    x_maxs = FT[];
    x_inis = FT[];
    Δ_inis = FT[];
    Δ_tols = FT[];
    if !isnothing(initial_guess)
        @assert length(initial_guess) == length(params);
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vcmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "b₆f" in params
            push!(x_mins, 0.01);
            push!(x_maxs, 5);
            iguess = findfirst(params .== "b₆f");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Γstar25" in params
            push!(x_mins, 1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Γstar25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Rd25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    else
        # loop through the data once to get the respiration and Γ_star limits
        rd_lim_min = 0.1;
        γ_lim_max = 40;
        for dfr in eachrow(df)
            # estimate the minimum respiration rate limit
            rd_lim_min = nanmax([rd_lim_min, -dfr.A_NET / temperature_correction(ps.trait.TD_R, dfr.T_LEAF)]);
            # estimate the maximum Γ_star limit
            γ_lim_max = nanmin([γ_lim_max, dfr.P_I / temperature_correction(ps.trait.TD_Γ, dfr.T_LEAF)]);
        end;
        # loop through the data once again to guess the Vcmax and b6f
        vcmax_guess = 5;
        b6f_guess = 0.1;
        ps.trait.r_d25 = rd_lim_min * 1.2;
        ps.trait.TD_Γ.VAL_REF = γ_lim_max * 0.8;
        for dfr in eachrow(df)
            photosynthesis!(config, ps, air, dfr.P_I, dfr.PPAR, dfr.T_LEAF);
            vcmax = (dfr.A_NET + ps.auxil.r_d) * (dfr.P_I + ps.auxil.k_m) / (dfr.P_I - ps.auxil.γ_star) / temperature_correction(ps.trait.TD_VCMAX, dfr.T_LEAF);
            vcmax_guess = nanmax([vcmax_guess, vcmax]);
            vcmax_guess = nanmin([vcmax_guess, 100]);
            j_1 = (dfr.A_NET + ps.auxil.r_d) * (4*dfr.P_I + 8*ps.auxil.γ_star) / (dfr.P_I - ps.auxil.γ_star) * ps.auxil.η;
            p_1 = dfr.PPAR * 0.5 * ps.auxil.ϕ_psi_max;
            v_q = p_1 * j_1 / max(1, p_1 - j_1);
            b6f_guess = nanmax([b6f_guess, v_q / ps.auxil.k_q]);
            b6f_guess = nanmin([b6f_guess, 1]);
        end;
        # set the initial guess
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, vcmax_guess * 2.0);
            push!(x_inis, vcmax_guess * 1.0);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "b₆f" in params
            push!(x_mins, 0.01);
            push!(x_maxs, b6f_guess * 2.0);
            push!(x_inis, b6f_guess * 1.0);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Γstar25" in params
            push!(x_mins, 1);
            push!(x_maxs, γ_lim_max * 1.0);
            push!(x_inis, γ_lim_max * 0.8);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Rd25" in params
            push!(x_mins, rd_lim_min * 1.0);
            push!(x_maxs, 5);
            push!(x_inis, rd_lim_min * 1.2);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    end;
    mthd = ReduceStepMethodND{FT}(x_mins = x_mins, x_maxs = x_maxs, x_inis = x_inis, Δ_inis = Δ_inis);
    stol = SolutionToleranceND{FT}(Δ_tols, 50);
    # func(x) = (rme = aci_rmse(config, ps, pst, air, df, x); @info "C3JB model" x rme; -rme);
    func(x) = -aci_rmse(config, ps, pst, air, df, params, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, params, sol);
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
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    x_mins = FT[];
    x_maxs = FT[];
    x_inis = FT[];
    Δ_inis = FT[];
    Δ_tols = FT[];
    if !isnothing(initial_guess)
        @assert length(initial_guess) == length(params);
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vcmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Rd25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    else
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            push!(x_inis, 50);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            push!(x_inis, 1);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    end;

    mthd = ReduceStepMethodND{FT}(x_mins = x_mins, x_maxs = x_maxs, x_inis = x_inis, Δ_inis = Δ_inis);
    stol = SolutionToleranceND{FT}(Δ_tols, 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, params, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, params, sol);
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
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    x_mins = FT[];
    x_maxs = FT[];
    x_inis = FT[];
    Δ_inis = FT[];
    Δ_tols = FT[];
    if !isnothing(initial_guess)
        @assert length(initial_guess) == length(params);
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vcmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Vpmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vpmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Rd25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    else
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            push!(x_inis, 50);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Vpmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            push!(x_inis, 50);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            push!(x_inis, 1);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    end;

    mthd = ReduceStepMethodND{FT}(x_mins = x_mins, x_maxs = x_maxs, x_inis = x_inis, Δ_inis = Δ_inis);
    stol = SolutionToleranceND{FT}(Δ_tols, 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, params, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, params, sol);
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
                params::Vector{String},
                initial_guess::Union{Nothing, Vector};
                min_count::Int = 9,
                rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve by removing outliers, given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `initial_guess` Initial guess of fitting parameters
- `min_count` Minimum number of data points to fit an A-Ci curve
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit_exclude_outliter(
            config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            params::Vector{String},
            initial_guess::Union{Nothing, Vector};
            min_count::Int = 9,
            rmse_threshold::Number = 2) where {FT}
    # remove outliers using thresholds when necessary
    df[!,"A_NET_BAK"] .= df.A_NET;
    last_rmse = 1000;
    last_sol = nothing;
    last_df = deepcopy(df);
    crnt_df = deepcopy(df);
    while true
        sol, best_rmse, aci = aci_fit(config, ps, air, crnt_df, params, initial_guess);
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
    best_rmse = aci_rmse(config, ps, air, df, params, last_sol);
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

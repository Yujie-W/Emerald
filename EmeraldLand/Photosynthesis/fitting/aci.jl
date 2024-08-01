#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain net photosynthetic rate
#
#######################################################################################################################################################################################################
"""

    aci_an(ps::LeafPhotosystem, air::AirLayer, p_i::Number, ppar::Number, t::Number)

Compute the net photosynthetic rate, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `ppar` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `t` Leaf temperature in `K`

"""
function aci_an(config::SPACConfiguration, ps::LeafPhotosystem, air::AirLayer, p_i::Number, ppar::Number, t::Number)
    photosynthesis!(config, ps, air, p_i, ppar, t);

    return ps.auxil.a_n
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain net photosynthetic rates of given data
#
#######################################################################################################################################################################################################
"""

    aci_curve(ps::LeafPhotosystem, air::AirLayer, pis::Vector, ppars::Vector, ts::Vector)
    aci_curve(ps::LeafPhotosystem, air::AirLayer, df::DataFrame)

Compute the net photosynthetic rates, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `pis` Internal CO₂ partial pressure in `Pa`
- `ppars` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `ts` Leaf temperature in `K`
- `df` DataFrame with columns `P_I`, `PPAR`, and `T_LEAF`

"""
function aci_curve end;

aci_curve(config::SPACConfiguration, ps::LeafPhotosystem, air::AirLayer, pis::Vector, ppars::Vector, ts::Vector) = aci_an.((config,), (ps,), (air,), pis, ppars, ts);

aci_curve(config::SPACConfiguration, ps::LeafPhotosystem, air::AirLayer, df::DataFrame) = aci_curve(config, ps, air, df.P_I, df.PPAR, df.T_LEAF);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain RMSE of A-Ci curve
#     2024-Jul-22: add support to C3CLM, C3FvCB, and C3VJP
#     2024-Jul-27: add option to turn on/off Rd fitting
#     2024-Jul-27: fit Γ_star as well for C3 models
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#
#######################################################################################################################################################################################################
"""

    aci_rmse(ps::LeafPhotosystem, air::AirLayer, df::DataFrame, params::Vector)

Compute the RMSE of A-Ci curve, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of parameters (Vcmax25, Vpmax25, b₆f, Jmax25, Rd25 depending on the photosynthesis model)

"""
function aci_rmse end;

aci_rmse(config::SPACConfiguration,
         ps::LeafPhotosystem,
         air::AirLayer,
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) = aci_rmse(config, ps, ps.trait, air, df, params; fit_rd = fit_rd, fitted_γ = fitted_γ);

aci_rmse(config::SPACConfiguration,
         ps::LeafPhotosystem,
         pst::Union{GeneralC3Trait, GeneralC4Trait},
         air::AirLayer,
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) = aci_rmse(config, ps, pst, pst.ACM, pst.AJM, pst.APM, air, df, params; fit_rd = fit_rd, fitted_γ = fitted_γ);

aci_rmse(config::SPACConfiguration,
         ps::LeafPhotosystem,
         pst::GeneralC3Trait,
         acm::AcMethodC3VcmaxPi,
         ajm::AjMethodC3JmaxPi,
         apm::ApMethodC3Vcmax,
         air::AirLayer,
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) = (
    pst.v_cmax25 = params[1];
    pst.j_max25 = params[2];
    if isnothing(fitted_γ)
        pst.TD_Γ.VAL_REF = params[3];
    else
        pst.TD_Γ.VAL_REF = fitted_γ;
    end;
    pst.r_d25 = fit_rd ? params[4] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);

aci_rmse(config::SPACConfiguration,
         ps::LeafPhotosystem,
         pst::GeneralC3Trait,
         acm::AcMethodC3VcmaxPi,
         ajm::AjMethodC3VqmaxPi,
         apm::ApMethodC3Vcmax,
         air::AirLayer,
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) = (
    pst.v_cmax25 = params[1];
    pst.b₆f = params[2];
    if isnothing(fitted_γ)
        pst.TD_Γ.VAL_REF = params[3];
    else
        pst.TD_Γ.VAL_REF = fitted_γ;
    end;
    pst.r_d25 = fit_rd ? params[4] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);

aci_rmse(config::SPACConfiguration,
         ps::LeafPhotosystem,
         pst::GeneralC4Trait,
         acm::AcMethodC4Vcmax,
         ajm::AjMethodC4JPSII,
         apm::ApMethodC4VcmaxPi,
         air::AirLayer,
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Nothing = nothing) = (
    pst.v_cmax25 = params[1];
    pst.r_d25 = fit_rd ? params[2] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);

aci_rmse(config::SPACConfiguration,
         ps::LeafPhotosystem,
         pst::GeneralC4Trait,
         acm::AcMethodC4Vcmax,
         ajm::AjMethodC4JPSII,
         apm::ApMethodC4VpmaxPi,
         air::AirLayer,
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Nothing = nothing) = (
    pst.v_cmax25 = params[1];
    pst.v_pmax25 = params[2];
    pst.r_d25 = fit_rd ? params[3] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);


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
#
#######################################################################################################################################################################################################
"""

    aci_fit(ps::LeafPhotosystem, air::AirLayer, df::DataFrame)

Fit the A-Ci curve, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`

"""
function aci_fit end;

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        air::AirLayer{FT},
        df::DataFrame;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing,
        initial_guess::Vector = [50, 2.1, 4, 1]) where {FT} = aci_fit(config, ps, ps.trait, air, df; fit_rd = fit_rd, fitted_γ = fitted_γ, initial_guess = initial_guess);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
        air::AirLayer,
        df::DataFrame;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing,
        initial_guess::Vector = [50, 2.1, 4, 1]) where {FT} = aci_fit(config, ps, pst, pst.ACM, pst.AJM, pst.APM, air, df; fit_rd = fit_rd, fitted_γ = fitted_γ, initial_guess = initial_guess);

aci_fit(config::SPACConfiguration{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC3Trait{FT},
        acm::AcMethodC3VcmaxPi,
        ajm::AjMethodC3JmaxPi,
        apm::ApMethodC3Vcmax,
        air::AirLayer,
        df::DataFrame;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing,
        initial_guess::Vector = [50, 100, 4, 1]) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 1, 1, 0.1] : [1, 1, 1],
        x_maxs = fit_rd ? [999, 999, 10, 10] : [999, 999, 10],
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
        air::AirLayer,
        df::DataFrame;
        fit_rd::Bool = false,
        fitted_γ::Union{Nothing, Number} = nothing,
        initial_guess::Vector = [50, 2.1, 4, 1]) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 0.01, 1, 0.1] : [1, 0.01, 1],
        x_maxs = fit_rd ? [999, 99, 10, 10] : [999, 99, 1],
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
        air::AirLayer,
        df::DataFrame;
        fit_rd::Bool = false,
        fitted_γ::Nothing = nothing,
        initial_guess::Vector = [50, 2.1, 4, 1]) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 0.1] : [1],
        x_maxs = fit_rd ? [999, 10] : [999],
        x_inis = fit_rd ? [100, 1] : [100],
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
        air::AirLayer,
        df::DataFrame;
        fit_rd::Bool = false,
        fitted_γ::Nothing = nothing,
        initial_guess::Vector = [50, 2.1, 4, 1]) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = fit_rd ? [1, 1, 0.1] : [1, 1],
        x_maxs = fit_rd ? [999, 999, 10] : [999, 999],
        x_inis = fit_rd ? [100, 100, 1] : [100, 100],
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
#
#######################################################################################################################################################################################################
"""

    aci_fit_exclude_outliter(ps::LeafPhotosystem{FT}, air::AirLayer{FT}, df::DataFrame; min_count::Int = 9, rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve by removing outliers, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `min_count` Minimum number of data points to fit an A-Ci curve
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit_exclude_outliter(
            config::SPACConfiguration{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame;
            fit_rd::Bool = false,
            fitted_γ::Union{Nothing, Number} = nothing,
            initial_guess::Vector = [50, 100, 4, 1],
            min_count::Int = 9,
            rmse_threshold::Number = 2) where {FT}
    # remove outliers using thresholds when necessary
    df[!,"A_NET_BAK"] .= df.A_NET;
    last_rmse = 1000;
    last_sol = nothing;
    last_df = deepcopy(df);
    crnt_df = deepcopy(df);
    while true
        sol, best_rmse, aci = aci_fit(config, ps, air, crnt_df; fit_rd = fit_rd, fitted_γ = fitted_γ, initial_guess = initial_guess);
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


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-20: add alias function to fit A-Ci curve
#     2024-Jul-23: add support to C3CytoInfApTrait
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#
#######################################################################################################################################################################################################
"""

    aci_fit!(df::DataFrame, model::String; min_count::Int = 9, remove_outlier::Bool = false, rmse_threshold::Number = 2)

Fit the A-Ci curve, given
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `model` Photosynthesis model string (C3Cyto, C3VJP, C4CLM, C4VJP)
- `min_count` Minimum number of data points to fit an A-Ci curve
- `remove_outlier` Remove outliers or not
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit!(
            config::SPACConfiguration,
            df::DataFrame,
            model::String;
            fit_rd::Bool = false,
            fitted_γ::Union{Nothing, Number} = nothing,
            initial_guess = nothing,
            min_count::Int = 9,
            remove_outlier::Bool = false,
            rmse_threshold::Number = 2)
    # first of all, make sure the DataFrame has the required columns
    @assert all([n in names(df) for n in ["P_I", "PPAR", "T_LEAF", "A_NET"]]) "The DataFrame should have columns P_I, PPAR, T_LEAF, and A_NET!";
    @assert nanmin(df.T_LEAF) > 253.15 "The leaf temperature should be in Kelvin!";

    # create a leaf photosystem based on the model string
    if model == "C3Cyto"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3CytoInfAp"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3JB"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        ps.trait.TD_ηC = ηCTDJohnson(Float64);
        ps.trait.TD_ηL = ηLTDJohnson(Float64);
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3JBInfAp"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        ps.trait.TD_ηC = ηCTDJohnson(Float64);
        ps.trait.TD_ηL = ηLTDJohnson(Float64);
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3VJP"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        new_guess = [50, 100, 4, 1];
    elseif model == "C3VJPInfAp"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        new_guess = [50, 100, 4, 1];
    elseif model == "C3CLM"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        ps.trait.COLIMIT_CJ = ColimitCJCLMC3(Float64);
        ps.trait.COLIMIT_IP = ColimitIPCLM(Float64);
        new_guess = [50, 100, 4, 1];
    elseif model == "C3CLMInfAp"
        ps = LeafPhotosystem{Float64}(trait = GeneralC3Trait{Float64}(), state = C3State{Float64}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        ps.trait.COLIMIT_CJ = ColimitCJCLMC3(Float64);
        ps.trait.COLIMIT_IP = ColimitIPCLM(Float64);
        new_guess = [50, 100, 4, 1];
    elseif model == "C4CLM"
        ps = LeafPhotosystem{Float64}(trait = GeneralC4Trait{Float64}(), state = C4State{Float64}());
        ps.trait.ACM = AcMethodC4Vcmax();
        ps.trait.AJM = AjMethodC4JPSII();
        ps.trait.APM = ApMethodC4VcmaxPi();
        new_guess = [50, 100, 4, 1];
    elseif model == "C4CLMSmooth"
        ps = LeafPhotosystem{Float64}(trait = GeneralC4Trait{Float64}(), state = C4State{Float64}());
        ps.trait.ACM = AcMethodC4Vcmax();
        ps.trait.AJM = AjMethodC4JPSII();
        ps.trait.APM = ApMethodC4VcmaxPi();
        ps.trait.COLIMIT_CJ = ColimitCJCLMC4(Float64);
        ps.trait.COLIMIT_IP = ColimitIPCLM(Float64);
        new_guess = [50, 100, 4, 1];
    elseif model == "C4VJP"
        ps = LeafPhotosystem{Float64}(trait = GeneralC4Trait{Float64}(), state = C4State{Float64}());
        ps.trait.ACM = AcMethodC4Vcmax();
        ps.trait.AJM = AjMethodC4JPSII();
        ps.trait.APM = ApMethodC4VpmaxPi();
        new_guess = [50, 100, 4, 1];
    else
        throw(ArgumentError("The model $(model) is not supported!"));
    end;

    # fit the A-Ci curve with or without removing outliers
    if remove_outlier
        return aci_fit_exclude_outliter(
                    config,
                    ps,
                    AirLayer{Float64}(),
                    df;
                    fit_rd = fit_rd,
                    fitted_γ = fitted_γ,
                    initial_guess = (isnothing(initial_guess) ? new_guess : initial_guess),
                    min_count = min_count,
                    rmse_threshold = rmse_threshold)
    else
        return aci_fit(
                    config,
                    ps,
                    AirLayer{Float64}(),
                    df;
                    fit_rd = fit_rd,
                    fitted_γ = fitted_γ,
                    initial_guess = (isnothing(initial_guess) ? new_guess : initial_guess))
    end;
end;






# TODO: move this function to EmeraldMath (not sure where yet)
function max_diff_index(y::Vector, p_pred::Vector)
    abs_diff = abs.(y .- p_pred);
    abs_diff[isnan.(abs_diff)] .= 0;

    return findmax(abs_diff)[2]
end;

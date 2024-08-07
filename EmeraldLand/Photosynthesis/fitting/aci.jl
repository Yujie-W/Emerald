#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-20: add alias function to fit A-Ci curve
#     2024-Jul-23: add support to C3CytoInfApTrait
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#     2024-Aug-06: move constructor to Namespace
#
#######################################################################################################################################################################################################
"""

    aci_fit!(config::SPACConfiguration{FT},
             df::DataFrame,
             model::String;
             fit_rd::Bool = false,
             fitted_γ::Union{Nothing, Number} = nothing,
             initial_guess::Union{Nothing, Vector} = nothing,
             min_count::Int = 9,
             remove_outlier::Bool = false,
             rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve, given
- `config` `SPACConfiguration` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `model` Photosynthesis model string (C3Cyto, C3CytoInfAp, C3JB, C3JBInfAp, C3VJP, C3VJPInfAp, C3CLM, C3CLMInfAp, C4CLM, C4CLMSmooth, C4VJP)
- `fit_rd` Fit Rd or not (default: false)
- `fitted_γ` Fitted Γ* value, will be prescribed if given (default: nothing)
- `initial_guess` Initial guess of fitting parameters
- `min_count` Minimum number of data points to fit an A-Ci curve
- `remove_outlier` Remove outliers or not
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit!(
            config::SPACConfiguration{FT},
            df::DataFrame,
            model::String;
            fit_rd::Bool = false,
            fitted_γ::Union{Nothing, Number} = nothing,
            initial_guess::Union{Nothing, Vector} = nothing,
            min_count::Int = 9,
            remove_outlier::Bool = false,
            rmse_threshold::Number = 2) where {FT}
    # first of all, make sure the DataFrame has the required columns
    @assert all([n in names(df) for n in ["P_I", "PPAR", "T_LEAF", "A_NET"]]) "The DataFrame should have columns P_I, PPAR, T_LEAF, and A_NET!";
    @assert nanmin(df.T_LEAF) > 253.15 "The leaf temperature should be in Kelvin!";

    # create a leaf photosystem based on the model string
    ps = LeafPhotosystem{FT}(model);
    new_guess =
        if model == "C3Cyto"
            [50, 2.1, 4, 1]             # Vcmax,  B6F, Gamma, Rd
        elseif model == "C3CytoInfAp"
            [50, 2.1, 4, 1]             # Vcmax,  B6F, Gamma, Rd
        elseif model == "C3JB"
            [50, 2.1, 4, 1]             # Vcmax,  B6F, Gamma, Rd
        elseif model == "C3JBInfAp"
            [50, 2.1, 4, 1]             # Vcmax,  B6F, Gamma, Rd
        elseif model == "C3VJP"
            [50, 100, 4, 1]             # Vcmax, Jmax, Gamma, Rd
        elseif model == "C3VJPInfAp"
            [50, 100, 4, 1]             # Vcmax, Jmax, Gamma, Rd
        elseif model == "C3CLM"
            [50, 100, 4, 1]             # Vcmax, Jmax, Gamma, Rd
        elseif model == "C3CLMInfAp"
            [50, 100, 4, 1]             # Vcmax, Jmax, Gamma, Rd
        elseif model == "C4CLM"
            [50, 1]                     # Vcmax, Rd
        elseif model == "C4CLMSmooth"
            [50, 1]                     # Vcmax, Rd
        elseif model == "C4VJP"
            [50, 50, 1]                 # Vcmax, Vpmax, Rd
        else
            error("The model $(model) is not supported!")
        end;

    # fit the A-Ci curve with or without removing outliers
    if remove_outlier
        return aci_fit_exclude_outliter(
                    config,
                    ps,
                    AirLayer{FT}(),
                    df,
                    (isnothing(initial_guess) ? new_guess : initial_guess);
                    fit_rd = fit_rd,
                    fitted_γ = fitted_γ,
                    min_count = min_count,
                    rmse_threshold = rmse_threshold)
    else
        return aci_fit(
                    config,
                    ps,
                    AirLayer{FT}(),
                    df,
                    (isnothing(initial_guess) ? new_guess : initial_guess);
                    fit_rd = fit_rd,
                    fitted_γ = fitted_γ)
    end;
end;

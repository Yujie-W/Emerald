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

    aci_fit!(config::SPACConfiguration{FT},
             df::DataFrame,
             model::String;
             fit_rd::Bool = false,
             fitted_γ::Union{Nothing, Number} = nothing,
             initial_guess = nothing,
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
            initial_guess = nothing,
            min_count::Int = 9,
            remove_outlier::Bool = false,
            rmse_threshold::Number = 2) where {FT}
    # first of all, make sure the DataFrame has the required columns
    @assert all([n in names(df) for n in ["P_I", "PPAR", "T_LEAF", "A_NET"]]) "The DataFrame should have columns P_I, PPAR, T_LEAF, and A_NET!";
    @assert nanmin(df.T_LEAF) > 253.15 "The leaf temperature should be in Kelvin!";

    # create a leaf photosystem based on the model string
    if model == "C3Cyto"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3CytoInfAp"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3JB"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        ps.trait.TD_ηC = ηCTDJohnson(FT);
        ps.trait.TD_ηL = ηLTDJohnson(FT);
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3JBInfAp"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3VqmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        ps.trait.TD_ηC = ηCTDJohnson(FT);
        ps.trait.TD_ηL = ηLTDJohnson(FT);
        new_guess = [50, 2.1, 4, 1];
    elseif model == "C3VJP"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        new_guess = [50, 100, 4, 1];
    elseif model == "C3VJPInfAp"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        new_guess = [50, 100, 4, 1];
    elseif model == "C3CLM"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Vcmax();
        ps.trait.COLIMIT_CJ = ColimitCJCLMC3(FT);
        ps.trait.COLIMIT_IP = ColimitIPCLM(FT);
        new_guess = [50, 100, 4, 1];
    elseif model == "C3CLMInfAp"
        ps = LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}());
        ps.trait.ACM = AcMethodC3VcmaxPi();
        ps.trait.AJM = AjMethodC3JmaxPi();
        ps.trait.APM = ApMethodC3Inf();
        ps.trait.COLIMIT_CJ = ColimitCJCLMC3(FT);
        ps.trait.COLIMIT_IP = ColimitIPCLM(FT);
        new_guess = [50, 100, 4, 1];
    elseif model == "C4CLM"
        ps = LeafPhotosystem{FT}(trait = GeneralC4Trait{FT}(), state = C4State{FT}());
        ps.trait.ACM = AcMethodC4Vcmax();
        ps.trait.AJM = AjMethodC4JPSII();
        ps.trait.APM = ApMethodC4VcmaxPi();
        new_guess = [50, 100, 4, 1];
    elseif model == "C4CLMSmooth"
        ps = LeafPhotosystem{FT}(trait = GeneralC4Trait{FT}(), state = C4State{FT}());
        ps.trait.ACM = AcMethodC4Vcmax();
        ps.trait.AJM = AjMethodC4JPSII();
        ps.trait.APM = ApMethodC4VcmaxPi();
        ps.trait.COLIMIT_CJ = ColimitCJCLMC4(FT);
        ps.trait.COLIMIT_IP = ColimitIPCLM(FT);
        new_guess = [50, 100, 4, 1];
    elseif model == "C4VJP"
        ps = LeafPhotosystem{FT}(trait = GeneralC4Trait{FT}(), state = C4State{FT}());
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
                    AirLayer{FT}(),
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
                    AirLayer{FT}(),
                    df;
                    fit_rd = fit_rd,
                    fitted_γ = fitted_γ,
                    initial_guess = (isnothing(initial_guess) ? new_guess : initial_guess))
    end;
end;

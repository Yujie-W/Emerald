#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-20: add alias function to fit A-Ci curve
#     2024-Jul-23: add support to C3CytoInfApTrait
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#     2024-Aug-06: move constructor to Namespace
#     2024-Oct-03: remove options for Rd and Γ_star fitting
#     2025-Jul-30: improved curve fitting algorithm
#
#######################################################################################################################################################################################################
"""

    aci_fit!(config::SPACConfiguration{FT},
             df::DataFrame,
             model::String,
             params::Vector{String};
             initial_guess::Union{Nothing, Vector} = nothing,
             min_count::Int = 9,
             remove_outlier::Bool = false,
             rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve, given
- `config` `SPACConfiguration` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `model` Photosynthesis model string (C3Cyto, C3CytoInfAp, C3JB, C3JBInfAp, C3VJP, C3VJPInfAp, C3CLM, C3CLMInfAp, C4CLM, C4CLMSmooth, C4VJP)
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `initial_guess` Initial guess of fitting parameters
- `min_count` Minimum number of data points to fit an A-Ci curve
- `remove_outlier` Remove outliers or not
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit!(
            config::SPACConfiguration{FT},
            df::DataFrame,
            model::String,
            params::Vector{String};
            initial_guess::Union{Nothing, Vector} = nothing,
            min_count::Int = 9,
            remove_outlier::Bool = false,
            rmse_threshold::Number = 2) where {FT}
    # first of all, make sure the DataFrame has the required columns
    @assert all([n in names(df) for n in ["P_I", "PPAR", "T_LEAF", "A_NET"]]) "The DataFrame should have columns P_I, PPAR, T_LEAF, and A_NET!";
    @assert nanmin(df.T_LEAF) > 253.15 "The leaf temperature should be in Kelvin!";

    # create a leaf photosystem based on the model string
    ps = LeafPhotosystem{FT}(model);

    # fit the A-Ci curve with or without removing outliers
    if remove_outlier
        return aci_fit_exclude_outliter(
                    config,
                    ps,
                    AirLayer{FT}(),
                    df,
                    params,
                    initial_guess;
                    min_count = min_count,
                    rmse_threshold = rmse_threshold)
    else
        return aci_fit(
                    config,
                    ps,
                    AirLayer{FT}(),
                    df,
                    params,
                    initial_guess)
    end;
end;

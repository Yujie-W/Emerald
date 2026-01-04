"""

    aci_fit!(config::SPACConfig{FT},
             df::DataFrame,
             c3c4::String,
             params::Vector{String};
             initial_guess::Union{Nothing, Vector} = nothing,
             min_count::Int = 9,
             remove_outlier::Bool = false,
             rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve, given
- `config` `SPACConfig` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `c3c4` C3 or C4 photosynthesis model
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `initial_guess` Initial guess of fitting parameters
- `min_count` Minimum number of data points to fit an A-Ci curve
- `remove_outlier` Remove outliers or not
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit!(
            config::SPACConfig{FT},
            df::DataFrame,
            c3c4::String,
            params::Vector{String};
            initial_guess::Union{Nothing, Vector} = nothing,
            min_count::Int = 9,
            remove_outlier::Bool = false,
            rmse_threshold::Number = 2) where {FT}
    # first of all, make sure the DataFrame has the required columns
    @assert all([n in names(df) for n in ["P_I", "PPAR", "T_LEAF", "A_NET"]]) "The DataFrame should have columns P_I, PPAR, T_LEAF, and A_NET!";
    @assert nanmin(df.T_LEAF) > 253.15 "The leaf temperature should be in Kelvin!";
    @assert c3c4 in ["C3", "C4"] "The c3c4 string should be either C3 or C4!";

    # create a leaf photosystem based on the c3c4 string
    cache = SPACCache{FT}(
                config.DIMENSIONS.DIM_AZI,
                config.DIMENSIONS.DIM_INCL,
                1,
                0,
                length(config.CONSTANTS.SPECTRA.Λ_SIF),
                length(config.CONSTANTS.SPECTRA.Λ_SIFE),
                length(config.CONSTANTS.SPECTRA.Λ));
    ps = LeafPhotosystem{FT}(c3c4);

    # fit the A-Ci curve with or without removing outliers
    return if remove_outlier
        aci_fit_exclude_outliter(
            config,
            cache,
            ps,
            AirLayer{FT}(),
            df,
            params,
            initial_guess;
            min_count = min_count,
            rmse_threshold = rmse_threshold)
    else
        aci_fit(
            config,
            cache,
            ps,
            AirLayer{FT}(),
            df,
            params,
            initial_guess)
    end;
end;

"""

    aci_fit(config::SPACConfig{FT},
            cache::SPACCache{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            df::DataFrame,
            params::Vector{String},
            initial_guess::Vector) where {FT}

Fit the A-Ci curve (will be abstractized based on the trait and methods embedded), given
- `config` `SPACConfig` struct
- `cache` `SPACCache` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `initial_guess` Initial guess of fitting parameters

"""
function aci_fit end;

aci_fit(config::SPACConfig{FT},
        ps::LeafPhotosystem{FT},
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = aci_fit(config, ps, ps.trait, air, df, params, initial_guess);

aci_fit(config::SPACConfig{FT},
        ps::LeafPhotosystem{FT},
        pst::C3Trait{FT},
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} =
    aci_fit(config, ps, pst, config.METHODS.C3_AC_METHOD, config.METHODS.C3_AJ_METHOD, config.METHODS.C3_AP_METHOD, air, df, params, initial_guess);

aci_fit(config::SPACConfig{FT},
        ps::LeafPhotosystem{FT},
        pst::C4Trait{FT},
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} =
    aci_fit(config, ps, pst, config.METHODS.C4_AC_METHOD, config.METHODS.C4_AJ_METHOD, config.METHODS.C4_AP_METHOD, air, df, params, initial_guess);


"""

    aci_fit_exclude_outliter(
                config::SPACConfig{FT},
                cache::SPACCache{FT},
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                df::DataFrame,
                params::Vector{String},
                initial_guess::Union{Nothing, Vector};
                min_count::Int = 9,
                rmse_threshold::Number = 2) where {FT}

Fit the A-Ci curve by removing outliers, given
- `config` `SPACConfig` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `initial_guess` Initial guess of fitting parameters
- `min_count` Minimum number of data points to fit an A-Ci curve
- `rmse_threshold` Threshold of RMSE to stop removing outliers

"""
function aci_fit_exclude_outliter(
            config::SPACConfig{FT},
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


max_diff_index(y::Vector, p_pred::Vector) = (
    abs_diff = abs.(y .- p_pred);
    abs_diff[isnan.(abs_diff)] .= 0;

    return findmax(abs_diff)[2]
);

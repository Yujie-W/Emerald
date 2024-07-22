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
function aci_an(ps::LeafPhotosystem, air::AirLayer, p_i::Number, ppar::Number, t::Number)
    photosynthesis!(ps, air, p_i, ppar, t);

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

aci_curve(ps::LeafPhotosystem, air::AirLayer, pis::Vector, ppars::Vector, ts::Vector) = aci_an.((ps,), (air,), pis, ppars, ts);

aci_curve(ps::LeafPhotosystem, air::AirLayer, df::DataFrame) = aci_curve(ps, air, df.P_I, df.PPAR, df.T_LEAF);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain RMSE of A-Ci curve
#     2024-Jul-22: add support to C3CLM, C3FvCB, and C3VJP
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

aci_rmse(ps::LeafPhotosystem, air::AirLayer, df::DataFrame, params::Vector) = aci_rmse(ps, ps.trait, air, df, params);

aci_rmse(ps::LeafPhotosystem, pst::C3CytoTrait, air::AirLayer, df::DataFrame, params::Vector) = (
    @assert length(params) == 3 "The number of parameters should be 3: Vcmax25, b₆f, Rd25!";
    pst.v_cmax25 = params[1];
    pst.b₆f = params[2];
    pst.r_d25 = params[3];

    return rmse(aci_curve(ps, air, df), df.A_NET)
);

aci_rmse(ps::LeafPhotosystem, pst::Union{C3CLMTrait, C3FvCBTrait, C3VJPTrait}, air::AirLayer, df::DataFrame, params::Vector) = (
    @assert length(params) == 3 "The number of parameters should be 3: Vcmax25, Jmax25, Rd25!";
    pst.v_cmax25 = params[1];
    pst.j_max25 = params[2];
    pst.r_d25 = params[3];

    return rmse(aci_curve(ps, air, df), df.A_NET)
);

aci_rmse(ps::LeafPhotosystem, pst::C4CLMTrait, air::AirLayer, df::DataFrame, params::Vector) = (
    @assert length(params) == 2 "The number of parameters should be 2: Vcmax25, Rd25!";
    pst.v_cmax25 = params[1];
    pst.r_d25 = params[2];

    return rmse(aci_curve(ps, air, df), df.A_NET)
);

aci_rmse(ps::LeafPhotosystem, pst::C4VJPTrait, air::AirLayer, df::DataFrame, params::Vector) = (
    @assert length(params) == 3 "The number of parameters should be 3: Vcmax25, Vpmax25, Rd25!";
    pst.v_cmax25 = params[1];
    pst.v_pmax25 = params[2];
    pst.r_d25 = params[3];

    return rmse(aci_curve(ps, air, df), df.A_NET)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to fit A-Ci curve
#     2024-Jul-22: add support to C3CLM, C3FvCB, and C3VJP
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

aci_fit(ps::LeafPhotosystem, air::AirLayer, df::DataFrame) = aci_fit(ps, ps.trait, air, df);

aci_fit(ps::LeafPhotosystem{FT}, pst::C3CytoTrait{FT}, air::AirLayer, df::DataFrame) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 0.01, 0.1],
        x_maxs = [300, 3, 10],
        x_inis = [50, 0.8, 1],
        Δ_inis = [10, 0.1, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.01, 0.01], 50);
    func(x) = -aci_rmse(ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(ps, pst, air, df, sol);
    aci = aci_curve(ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(ps::LeafPhotosystem{FT}, pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}}, air::AirLayer, df::DataFrame) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 1, 0.1],
        x_maxs = [300, 600, 10],
        x_inis = [50, 100, 1],
        Δ_inis = [10, 10, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.1, 0.01], 50);
    func(x) = -aci_rmse(ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(ps, pst, air, df, sol);
    aci = aci_curve(ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(ps::LeafPhotosystem{FT}, pst::C4CLMTrait{FT}, air::AirLayer, df::DataFrame) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 0.1],
        x_maxs = [200, 10],
        x_inis = [100, 1],
        Δ_inis = [10, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.01], 50);
    func(x) = -aci_rmse(ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(ps, pst, air, df, sol);
    aci = aci_curve(ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(ps::LeafPhotosystem{FT}, pst::C4VJPTrait{FT}, air::AirLayer, df::DataFrame) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 1, 0.1],
        x_maxs = [200, 200, 10],
        x_inis = [100, 100, 1],
        Δ_inis = [10, 10, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.1, 0.01], 50);
    func(x) = -aci_rmse(ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(ps, pst, air, df, sol);
    aci = aci_curve(ps, air, df);

    return sol, best_rmse, aci
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-20: add functions to fit A-Ci curve with removing outliers
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
function aci_fit_exclude_outliter(ps::LeafPhotosystem{FT}, air::AirLayer{FT}, df::DataFrame; min_count::Int = 9, rmse_threshold::Number = 2) where {FT}
    # remove outliers using thresholds when necessary
    df[!,"A_NET_BAK"] .= df.A_NET;
    last_rmse = 1000;
    last_sol = nothing;
    last_df = deepcopy(df);
    crnt_df = deepcopy(df);
    while true
        sol, best_rmse, aci = aci_fit(ps, air, crnt_df);
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
    best_rmse = aci_rmse(ps, air, df, last_sol);
    aci = aci_curve(ps, air, df);

    return last_sol, best_rmse, aci
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-20: add alias function to fit A-Ci curve
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
function aci_fit!(df::DataFrame, model::String; min_count::Int = 9, remove_outlier::Bool = false, rmse_threshold::Number = 2)
    # first of all, make sure the DataFrame has the required columns
    @assert all([n in names(df) for n in ["P_I", "PPAR", "T_LEAF", "A_NET"]]) "The DataFrame should have columns P_I, PPAR, T_LEAF, and A_NET!";
    @assert nanmin(df.T_LEAF) > 253.15 "The leaf temperature should be in Kelvin!";

    # create a leaf photosystem based on the model string
    ps = if model == "C3Cyto"
        LeafPhotosystem{Float64}(trait = C3CytoTrait{Float64}(), state = C3CytoState{Float64}())
    elseif model == "C3CLM"
        LeafPhotosystem{Float64}(trait = C3CLMTrait{Float64}(), state = C3VJPState{Float64}())
    elseif model == "C3FvCB"
        LeafPhotosystem{Float64}(trait = C3FvCBTrait{Float64}(), state = C3VJPState{Float64}())
    elseif model == "C3VJP"
        LeafPhotosystem{Float64}(trait = C3VJPTrait{Float64}(), state = C3VJPState{Float64}())
    elseif model == "C4CLM"
        LeafPhotosystem{Float64}(trait = C4CLMTrait{Float64}(), state = C4VJPState{Float64}())
    elseif model == "C4VJP"
        LeafPhotosystem{Float64}(trait = C4VJPTrait{Float64}(), state = C4VJPState{Float64}())
    else
        throw(ArgumentError("The model should be one of C3Cyto, C3VJP, C4CLM, and C4VJP!"))
    end;

    # fit the A-Ci curve with or without removing outliers
    if remove_outlier
        return aci_fit_exclude_outliter(ps, AirLayer{Float64}(), df, min_count = min_count, rmse_threshold = rmse_threshold)
    else
        return aci_fit(ps, AirLayer{Float64}(), df)
    end;
end;






# TODO: move this function to EmeraldMath (not sure where yet)
function max_diff_index(y::Vector, p_pred::Vector)
    abs_diff = abs.(y .- p_pred);
    abs_diff[isnan.(abs_diff)] .= 0;

    return findmax(abs_diff)[2]
end;

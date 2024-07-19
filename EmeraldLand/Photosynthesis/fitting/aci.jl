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

aci_curve(ps::LeafPhotosystem, air::AirLayer, df::DataFrame) = aci_an.(ps, air, df.P_I, df.PPAR, df.T_LEAF);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain RMSE of A-Ci curve
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

aci_rmse(ps::LeafPhotosystem, pst::C3VJPTrait, air::AirLayer, df::DataFrame, params::Vector) = (
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

aci_fit(ps::LeafPhotosystem, air::AirLayer, df::DataFrame) = nothing;

aci_fit(ps::LeafPhotosystem{FT}, pst::C3CytoTrait{FT}, air::AirLayer, df::DataFrame) where {FT} = (
    mthd = ReduceStepMethodND{FT}(
        x_mins = [1, 0.01, 0.1],
        x_maxs = [300, 3, 10],
        x_inis = [50, 0.8, 1],
        Δ_inis = [10, 0.1, 1],
    );
    stol = SolutionToleranceND{FT}([0.1, 0.1, 0.01], 50);
    func(x) = -aci_rmse(ps, pst, air, df, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(ps, pst, air, df, sol);
    aci = aci_curve(ps, air, df);

    return sol, best_rmse, aci
);

aci_fit(ps::LeafPhotosystem{FT}, pst::C3VJPTrait{FT}, air::AirLayer, df::DataFrame) where {FT} = (
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

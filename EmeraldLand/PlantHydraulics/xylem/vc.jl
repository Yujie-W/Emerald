# This file contains the functions to compute the vulnerability curve related variables

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Feb-01: add method for LogisticVC
#     2022-Feb-01: add method for PowerVC
#     2022-Feb-01: add method for WeibullVC
#     2022-Feb-01: add method for ComplexVC
#     2022-Oct-21: make sure relative K does not exceed 1
#     2024-Feb-29: set the minimum value of relative K to eps(FT)
#
#######################################################################################################################################################################################################
"""

    relative_xylem_k(vc::ComplexVC{FT}, p_25::FT) where {FT}
    relative_xylem_k(vc::LogisticVC{FT}, p_25::FT) where {FT}
    relative_xylem_k(vc::PowerVC{FT}, p_25::FT) where {FT}
    relative_xylem_k(vc::WeibullVC{FT}, p_25::FT) where {FT}

Return the relative xylem hydraulic conductance, given
- `vc` `ComplexVC`, `LogisticVC`, `PowerVC`, or `WeibullVC` type vulnerability curve
- `p_25` Equivalent xylem water pressure at 298.15 K in `[MPa]` (surface tension correction made already)

"""
function relative_xylem_k end;

relative_xylem_k(vc::ComplexVC{FT}, p_25::FT) where {FT} = (
    @assert sum(vc.fs) ≈ 1 "Probabilities of VCs must sum up to 1!";
    @assert length(vc.fs) == length(vc.vcs) "Lengths of VC curves and probabilities must equal!";

    kr::FT = 0;
    for i in eachindex(vc.fs)
        kr += relative_xylem_k(vc.vcs[i], p_25) * vc.fs[i];
    end;

    return kr
);

relative_xylem_k(vc::LogisticVC{FT}, p_25::FT) where {FT} = (
    if p_25 >= 0
        return FT(1)
    end;

    return min(1, max(eps(FT), (1 - 1/(1 + vc.A * exp(vc.B * p_25))) * (vc.A + 1) / vc.A))
);

relative_xylem_k(vc::PowerVC{FT}, p_25::FT) where {FT} = (
    if p_25 >= 0
        return FT(1)
    end;

    return min(1, max(eps(FT), 1 / (1 + vc.A * (-p_25) ^ vc.B)))
);

relative_xylem_k(vc::WeibullVC{FT}, p_25::FT) where {FT} = (
    if p_25 >= 0
        return FT(1)
    end;

    return min(1, max(eps(FT), exp(-1 * (-p_25 / vc.B) ^ vc.C)))
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Feb-02: add method for LogisticVC
#     2022-Feb-02: add method for PowerVC
#     2022-Feb-02: add method for WeibullVC
#     2022-Feb-02: add method for ComplexVC
#     2022-Feb-02: add a reference kr for more customized calculations
#     2023-Sep-22: use solver to find the target pressure for ComplexVC
# Bug fixes
#     2023-Mar-02: fix an issue with Weibull function critical pressure
#
#######################################################################################################################################################################################################
"""

    xylem_pressure(vc::ComplexVC{FT}, kr::FT) where {FT}
    xylem_pressure(vc::LogisticVC{FT}, kr::FT) where {FT}
    xylem_pressure(vc::PowerVC{FT}, kr::FT) where {FT}
    xylem_pressure(vc::WeibullVC{FT}, kr::FT) where {FT}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ComplexVC`, `LogisticVC`, `PowerVC`, or `WeibullVC` type struct
- `kr` Reference conductance

"""
function xylem_pressure end;

xylem_pressure(vc::LogisticVC{FT}, kr::FT) where {FT} = log(kr / (vc.A + 1 - kr * vc.A)) / vc.B;

xylem_pressure(vc::PowerVC{FT}, kr::FT) where {FT} = -1 * ((1 - kr) / (kr * vc.A)) ^ (1 / vc.B);

xylem_pressure(vc::WeibullVC{FT}, kr::FT) where {FT} = -1 * (-1 * log(kr)) ^ (1 / vc.C) * vc.B;

xylem_pressure(vc::ComplexVC{FT}, kr::FT) where {FT} = (
    # use the p from each curve to define the min and max of target pressure
    p_min::FT = 0;
    p_max::FT = -Inf;
    for sub in vc.vcs
        p = xylem_pressure(sub, kr);
        if p > p_max
            p_max = p;
        end;
        if p < p_min
            p_min = p;
        end;
    end;

    # use the Solver function to find the target pressure
    mthd = NewtonBisectionMethod{FT}(x_min=p_min, x_max=p_max, x_ini=(p_min+p_max)/2);
    stol = SolutionTolerance{FT}(eps(FT)*100, 50);
    @inline f(p) = relative_xylem_k(vc, p) - kr;

    return find_zero(f, mthd, stol)
);

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Feb-02: migrate function to version v0.3
#     2022-Feb-02: rename the function to critical_pressure
#     2022-Feb-02: add method for LogisticVC
#     2022-Feb-02: add method for PowerVC
#     2022-Feb-02: add method for WeibullVC
#     2022-Feb-02: add method for ComplexVC
#     2022-Feb-02: add a reference kr for more customized calculations
#     2022-May-25: iterate through VCS rather than its indices for ComplexVC
#     2023-Aug-27: add nan check
#     2023-Sep-11: remove the default value of kr
#     2023-Sep-11: renmae function to xylem_pressure to be general
# Bug fixes
#     2023-Mar-02: fix an issue with Weibull function critical pressure
#
#######################################################################################################################################################################################################
"""

    xylem_pressure(vc::ComplexVC{FT}, kr::FT) where {FT<:AbstractFloat}
    xylem_pressure(vc::LogisticVC{FT}, kr::FT) where {FT<:AbstractFloat}
    xylem_pressure(vc::PowerVC{FT}, kr::FT) where {FT<:AbstractFloat}
    xylem_pressure(vc::WeibullVC{FT}, kr::FT) where {FT<:AbstractFloat}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ComplexVC`, `LogisticVC`, `PowerVC`, or `WeibullVC` type struct
- `kr` Reference conductance

"""
function xylem_pressure end

xylem_pressure(vc::LogisticVC{FT}, kr::FT) where {FT<:AbstractFloat} = (
    _p = log(kr / (vc.A + 1 - kr * vc.A)) / vc.B;

    return isnan(_p) ? error("NaN detected when computing critial pressure!") : _p
);

xylem_pressure(vc::PowerVC{FT}, kr::FT) where {FT<:AbstractFloat} = (
    _p = -1 * ((1 - kr) / (kr * vc.A)) ^ (1 / vc.B);

    return isnan(_p) ? error("NaN detected when computing critial pressure!") : _p
);

xylem_pressure(vc::WeibullVC{FT}, kr::FT) where {FT<:AbstractFloat} = (
    _p = -1 * (-1 * log(kr)) ^ (1 / vc.C) * vc.B;

    return isnan(_p) ? error("NaN detected when computing critial pressure!") : _p
);

xylem_pressure(vc::ComplexVC{FT}, kr::FT) where {FT<:AbstractFloat} = (
    (; VCS) = vc;

    _p_crit::FT = 0;
    for _vc in VCS
        _p_crit = min(_p_crit, xylem_pressure(_vc, kr));
    end;

    return isnan(_p_crit) ? error("NaN detected when computing critial pressure!") : _p_crit
);

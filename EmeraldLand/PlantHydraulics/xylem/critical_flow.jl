# This file contains functions to compute leaf critical flow

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add function to compute xylem end; pressure at steady state (no water exchange through the capacitor along the xylem)
#     2024-Feb-28: add LAI <= 0 control
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(xylem::XylemHydraulics{FT}, flow::FT, t::FT) where {FT}

Return the xylem pressure at the end; of xylem, given
- `xylem` `XylemHydraulics` type struct
- `flow` flow rate through the xylem `[mol s⁻¹]`
- `t` Xylem water temperature `[K]`

"""
function xylem_end_pressure(xylem::XylemHydraulics{FT}, flow::FT, t::FT) where {FT}
    if xylem.trait.area <= 0
        return FT(0)
    end;

    # run the pressure profile calculation only if xylem area > 0
    k_max = xylem.trait.area * xylem.trait.k_max / xylem.trait.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(xylem.auxil.k_history);
    p = xylem.auxil.pressure[1];
    for i in 1:N
        p_mem = xylem.state.p_history[i];
        k_mem = xylem.auxil.k_history[i];

        p₂₅ = p / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(xylem.trait.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        # flow rate is the mean of that at two planes (i and i+1)
        p -= flow / k + ρg_MPa(FT) * xylem.trait.Δh / N;
    end;

    return p
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add function to compute xylem critical pressure at steady state (no water exchange through the capacitor along the xylem)
#     2024-Feb-28: add LAI <= 0 control
#     2024_Jul-24: use spac cache
#
#######################################################################################################################################################################################################
"""

    critical_flow(config::SPACConfiguration{FT}, xylem::XylemHydraulics{FT}, cache::SPACCache{FT}, t::FT, ini::FT = FT(0.5)) where {FT}

Return the critical flow rate that triggers a given amount of loss of hydraulic conductance, given
- `config` `SPACConfiguration` type struct
- `xylem` `XylemHydraulics` type struct
- `cache` `SPACCache` type struct
- `t` Xylem water temperature `[K]`
- `ini` Initial guess

"""
function critical_flow(config::SPACConfiguration{FT}, xylem::XylemHydraulics{FT}, cache::SPACCache{FT}, t::FT, ini::FT = FT(0.5)) where {FT}
    if xylem.trait.area <= 0
        return FT(0)
    end;

    # run the pressure profile calculation only if xylem area > 0
    (; KR_THRESHOLD) = config;

    # compute the misc variables
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);
    p_crt = xylem_pressure(xylem.trait.vc, KR_THRESHOLD) * f_st;

    # add a judgement to make sure p_ups is higher than p_crt
    if xylem.auxil.pressure[1] < p_crt
        return eps(FT)
    end;

    # set up method to calculate critical flow
    fh = (xylem.auxil.pressure[1] - p_crt) * xylem.trait.k_max * xylem.trait.area / f_vis;
    fl = FT(0);
    ms = cache.solver_nb;
    ms.x_min = fl;
    ms.x_max = fh;
    ms.x_ini = min((fh + fl) / 2, ini);
    st = cache.stol_nb;

    # define the target function
    @inline f(x) = xylem_end_pressure(xylem, x, t) - p_crt;

    # find the solution
    solut = find_zero(f, ms, st);

    # warning if the solution is NaN
    if isnan(solut)
        # @warn "E_crit is NaN, please check the settings..." hs.p_ups;
        solut = eps(FT);
    end;

    return solut
end;

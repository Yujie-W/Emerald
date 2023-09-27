#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-01: migrate the function
#     2022-Jun-01: add method for LeafHydraulics
#     2022-Oct-21: add a p_ups if statement
#     2023-Sep-11: add config to the variable list
#     2023-Sep-11: rename function to critical_flow
#
#######################################################################################################################################################################################################
"""
This function returns the critical flow rate that triggers a given amount of loss of hydraulic conductance for
- Leaf hydraulic system
- Mono element SPAC system

"""
function critical_flow end


"""

    critical_flow(config::SPACConfiguration{FT}, hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5)) where {FT}

Return the critical flow rate that triggers a given amount of loss of conductance, given
- `config` `SPACConfiguration` type struct
- `hs` `LeafHydraulics` type struct
- `T` Liquid temperature
- `ini` Initial guess

"""
critical_flow(config::SPACConfiguration{FT}, hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5)) where {FT} = (
    (; KR_THRESHOLD) = config;
    (; K_SLA, VC) = hs;

    # compute the misc variables
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_crt = xylem_pressure(VC, KR_THRESHOLD) * _f_st;

    # add a judgement to make sure p_ups is higher than _p_crt
    if hs.p_ups < _p_crt
        return eps(FT)
    end;

    # set up method to calculate critical flow
    _fh = (hs.p_ups - _p_crt) * K_SLA / _f_vis;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _st = SolutionTolerance{FT}(eps(FT)*100, 50);

    # define the target function
    @inline f(x) = xylem_end_pressure(hs, x, T) - _p_crt;

    # find the solution
    _solut = find_zero(f, _ms, _st);

    # warning if the solution is NaN
    if isnan(_solut)
        # @warn "E_crit is NaN, please check the settings..." hs.p_ups;
        _solut = eps(FT);
    end;

    return _solut
);

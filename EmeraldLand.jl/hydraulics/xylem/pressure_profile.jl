# This file contains the functions and methods used to update the pressure profile within the xylem

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function to update the xylem pressure profile
#
# #######################################################################################################################################################################################################
"""

    xylem_pressure_profile!(xylem::XylemHydraulics{FT}, t::FT) where {FT}

Update the xylem pressure profile, given
- `xylem` `XylemHydraulics` type struct
- `t` Xylem water temperature `[K]`

"""
function xylem_pressure_profile! end;

xylem_pressure_profile!(xylem::XylemHydraulics{FT}, t::FT) where {FT} = xylem_pressure_profile!(xylem.state, xylem.auxil, t);

xylem_pressure_profile!(x_state::XylemHydraulicsState{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, t::FT) where {FT} = (
    k_max = x_state.area * x_state.k_max / x_state.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(x_aux.k_history);
    for i in 1:N
        p_mem = x_state.p_history[i];
        k_mem = x_aux.k_history[i];

        p₂₅ = x_aux.pressure[i] / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(x_state.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        # flow rate is the mean of that at two planes (i and i+1)
        x_aux.pressure[i+1] -= (x_aux.flow[i] + x_aux.flow[i+1]) / 2 / k + ρg_MPa(FT) * x_state.Δh / N;
    end;

    return nothing
);

xylem_pressure_profile!(x_state::XylemHydraulicsState{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, t::FT) where {FT} = (
    k_max = x_state.area * x_state.k_max / x_state.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(x_aux.k_history);
    for i in 1:N
        p_mem = x_state.p_history[i];
        k_mem = x_aux.k_history[i];

        p₂₅ = x_aux.pressure[i] / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(x_state.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        x_aux.pressure[i+1] -= x_aux.flow / k + ρg_MPa(FT) * x_state.Δh / N;
    end;

    return nothing
);

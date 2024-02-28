# This file contains the functions and methods used to update the pressure profile within the xylem

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function to update the xylem pressure profile
#     2023-Sep-30: add method to update the xylem pressure profile in the reverse order
#     2024-Feb-28: add LAI <= 0 control
#
#######################################################################################################################################################################################################
"""

    xylem_pressure_profile!(xylem::XylemHydraulics{FT}, t::FT) where {FT}
    xylem_pressure_profile!(xylem::XylemHydraulics{FT}, t::FT, rev::Bool) where {FT}

Update the xylem pressure profile, given
- `xylem` `XylemHydraulics` type struct
- `t` Xylem water temperature `[K]`
- `rev` If `rev` is given (true or false), update the pressure profile in the reverse order from N to 1

"""
function xylem_pressure_profile! end;

xylem_pressure_profile!(xylem::XylemHydraulics{FT}, t::FT) where {FT} = xylem_pressure_profile!(xylem.trait, xylem.state, xylem.auxil, t);

xylem_pressure_profile!(x_trait::XylemHydraulicsTrait{FT}, x_state::XylemHydraulicsState{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, t::FT) where {FT} = (
    if x_trait.area <= 0
        return nothing
    end;

    # update the pressure profile calculation only if xylem area > 0
    k_max = x_trait.area * x_trait.k_max / x_trait.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(x_aux.k_history);
    for i in 1:N
        p_mem = x_state.p_history[i];
        k_mem = x_aux.k_history[i];

        p₂₅ = x_aux.pressure[i] / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(x_trait.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        # flow rate is the mean of that at two planes (i and i+1)
        x_aux.pressure[i+1] = x_aux.pressure[i] - (x_aux.flow[i] + x_aux.flow[i+1]) / 2 / k - ρg_MPa(FT) * x_trait.Δh / N;
    end;

    return nothing
);

xylem_pressure_profile!(x_trait::XylemHydraulicsTrait{FT}, x_state::XylemHydraulicsState{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, t::FT) where {FT} = (
    if x_trait.area <= 0
        return nothing
    end;

    # update the pressure profile calculation only if xylem area > 0
    k_max = x_trait.area * x_trait.k_max / x_trait.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(x_aux.k_history);
    for i in 1:N
        p_mem = x_state.p_history[i];
        k_mem = x_aux.k_history[i];

        p₂₅ = x_aux.pressure[i] / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(x_trait.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        x_aux.pressure[i+1] = x_aux.pressure[i] - x_aux.flow / k - ρg_MPa(FT) * x_trait.Δh / N;
    end;

    return nothing
);

xylem_pressure_profile!(xylem::XylemHydraulics{FT}, t::FT, rev::Bool) where {FT} = xylem_pressure_profile!(xylem.trait, xylem.state, xylem.auxil, t, rev);

xylem_pressure_profile!(x_trait::XylemHydraulicsTrait{FT}, x_state::XylemHydraulicsState{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, t::FT, ::Bool) where {FT} = (
    if x_trait.area <= 0
        return nothing
    end;

    # update the pressure profile calculation only if xylem area > 0
    k_max = x_trait.area * x_trait.k_max / x_trait.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(x_aux.k_history);
    for i in N:-1:1
        p_mem = x_state.p_history[i];
        k_mem = x_aux.k_history[i];

        p₂₅ = x_aux.pressure[i+1] / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(x_trait.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        # flow rate is the mean of that at two planes (i and i+1)
        x_aux.pressure[i] = x_aux.pressure[i+1] + (x_aux.flow[i] + x_aux.flow[i+1]) / 2 / k + ρg_MPa(FT) * x_trait.Δh / N;
    end;

    return nothing
);

xylem_pressure_profile!(x_trait::XylemHydraulicsTrait{FT}, x_state::XylemHydraulicsState{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, t::FT, ::Bool) where {FT} = (
    if x_trait.area <= 0
        return nothing
    end;

    # update the pressure profile calculation only if xylem area > 0
    k_max = x_trait.area * x_trait.k_max / x_trait.l;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);

    N = length(x_aux.k_history);
    for i in N:-1:1
        p_mem = x_state.p_history[i];
        k_mem = x_aux.k_history[i];

        p₂₅ = x_aux.pressure[i+1] / f_st;
        if p₂₅ < p_mem
            k = relative_xylem_k(x_trait.vc, p₂₅) / f_vis * k_max * N;
        else
            k = k_mem / f_vis * k_max * N;
        end;

        x_aux.pressure[i] = x_aux.pressure[i+1] + x_aux.flow / k + ρg_MPa(FT) * x_trait.Δh / N;
    end;

    return nothing
);

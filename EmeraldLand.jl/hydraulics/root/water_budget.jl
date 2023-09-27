# This file contains function to update the water budget of the root

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#
#######################################################################################################################################################################################################
"""

    root_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}

Set the flow profile of each root, given
- `spac` `MultiLayerSPAC` type struct
- `Δt` time step

"""
function root_water_budget! end

root_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT} = (
    for root in spac.ROOTS
        root_water_budget!(root.NS, root.NS.xylem.auxil, Δt);
    end;

    return nothing
);

root_water_budget!(root::Root{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, Δt::FT) where {FT} = (
    f_vis = relative_viscosity(root.energy.auxil.t);

    # make sure the buffer rate does not drain or overflow the capacictance
    N = root.xylem.state.v_storage;
    for i in 1:N
        if root.xylem.auxil.flow_buffer[i] > 0 && root.xylem.state.v_storage[i] <= root.xylem.auxil.flow_buffer[i] * Δt
            @warn "The capacitance buffer is drained, use only half of the remaining water in the buffer!";
            root.xylem.auxil.flow_buffer[i] = root.xylem.state.v_storage[i] / 2 / Δt;
        end;
    end;

    # update the integrators of the flow (do not use flow_out here as it may be higher than the flow_in + sum_flow_buffer)
    root.∫∂w∂t_in += flow_in(x_aux) * Δt;
    root.∫∂w∂t_out += (flow_in(x_aux) + sum(root.xylem.auxil.flow_buffer)) * Δt;

    # update storage and the tissue pressure (p_storage)
    v_max_i = root.xylem.state.v_max * root.xylem.state.area * root.xylem.state.l / N;
    for i in 1:N
        root.xylem.state.v_storage[i] -= root.xylem.auxil.flow_buffer[i] * Δt;
        root.xylem.state.p_storage[i] = capacitance_pressure(root.xylem.state.pv, root.xylem.state.v_storage[i] / v_max_i, root.energy.auxil.t);
        root.xylem.auxil.flow_buffer[i] = (root.xylem.state.p_storage[i] - (x_aux.pressure[i] + x_aux.pressure[i+1]) / 2) * root.xylem.state.pv.k_refill / f_vis * root.xylem.state.v_storage[i];
    end;

    return nothing
);

root_water_budget!(root::Root{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, Δt::FT) where {FT} = (
    # update the integrators of the flow
    root.∫∂w∂t_in += flow_in(x_aux) * Δt;
    root.∫∂w∂t_out += flow_out(x_aux) * Δt;

    return nothing
);

# This file contains function to update the water budget of the root

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#
#######################################################################################################################################################################################################
"""

    xylem_water_budget!(organ::Union{Root{FT}, Stem{FT}}, xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, t::FT, δt::FT) where {FT}
    xylem_water_budget!(organ::Union{Root{FT}, Stem{FT}}, xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, t::FT, δt::FT) where {FT}


Set the flow profile of the root or stem xylem, given
- `organ` `Root` or `Stem` type struct
- `xylem` `XylemHydraulics` type struct
- `x_aux` `XylemHydraulicsAuxilNSS` or `XylemHydraulicsAuxilSS` type struct
- `t` temperature `[K]`
- `δt` time step

"""
function xylem_water_budget! end

xylem_water_budget!(organ::Union{Root2{FT}, Stem2{FT}}, xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, t::FT, δt::FT) where {FT} = (
    f_vis = relative_viscosity(t);

    # make sure the buffer rate does not drain or overflow the capacictance
    N = length(xylem.state.v_storage);
    for i in 1:N
        if xylem.auxil.flow_buffer[i] > 0 && xylem.state.v_storage[i] <= xylem.auxil.flow_buffer[i] * δt
            @warn "The capacitance buffer is drained, use only half of the remaining water in the buffer!";
            xylem.auxil.flow_buffer[i] = xylem.state.v_storage[i] / 2 / δt;
        end;
    end;

    # update the integrators of the flow (do not use flow_out here as it may be higher than the flow_in + sum_flow_buffer)
    organ.∫∂w∂t_in += flow_in(xylem) * δt;
    organ.∫∂w∂t_out += (flow_in(xylem) + sum(xylem.auxil.flow_buffer)) * δt;

    # update storage and the tissue pressure (p_storage)
    v_max_i = xylem.state.v_max * xylem.state.area * xylem.state.l / N;
    for i in 1:N
        xylem.state.v_storage[i] -= xylem.auxil.flow_buffer[i] * δt;
        xylem.auxil.p_storage[i] = capacitance_pressure(xylem.state.pv, xylem.state.v_storage[i] / v_max_i, t);
        xylem.auxil.flow_buffer[i] = (xylem.auxil.p_storage[i] - (x_aux.pressure[i] + x_aux.pressure[i+1]) / 2) * xylem.state.pv.k_refill / f_vis * xylem.state.v_storage[i];
    end;

    return nothing
);

xylem_water_budget!(organ::Union{Root2{FT}, Stem2{FT}}, xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, t::FT, δt::FT) where {FT} = (
    # update the integrators of the flow
    organ.∫∂w∂t_in += flow_in(xylem) * δt;
    organ.∫∂w∂t_out += flow_out(xylem) * δt;

    return nothing
);

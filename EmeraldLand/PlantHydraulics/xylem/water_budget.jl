# This file contains function to update the water budget of the root

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add root_water_budget! function
#     2024-Feb-28: add LAI <= 0 control
#     2024-Sep-03: use state.asap to check the xylem status (<= 0 means the xylem is dead)
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
function xylem_water_budget! end;

xylem_water_budget!(xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, δt::FT) where {FT} = (
    if xylem.state.asap <= 0
        return nothing
    end;

    # run the water budget for the xylem only if xylem area > 0
    # make sure the buffer rate does not drain or overflow the capacictance
    N = length(xylem.state.v_storage);
    for i in 1:N
        if xylem.auxil.flow_buffer[i] > 0 && xylem.state.v_storage[i] <= xylem.auxil.flow_buffer[i] * δt
            @warn "The capacitance buffer is drained, use only half of the remaining water in the buffer!";
            xylem.auxil.flow_buffer[i] = xylem.state.v_storage[i] / 2 / δt;
        end;
        xylem.state.v_storage[i] -= xylem.auxil.flow_buffer[i] * δt;
    end;

    return nothing
);

xylem_water_budget!(xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, δt::FT) where {FT} = nothing;

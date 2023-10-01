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

xylem_water_budget!(xylem::XylemHydraulics{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, δt::FT) where {FT} = (
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

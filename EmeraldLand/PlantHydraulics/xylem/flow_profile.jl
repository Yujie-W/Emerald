# This function is used to read and set up the flow rate profile in a xylem element

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-22: add function to read the flow rate into the xylem
#
#######################################################################################################################################################################################################
"""

    flow_in(organ)

Return the flow rate into the organ

"""
function flow_in end;

flow_in(xylem::XylemHydraulics{FT}) where {FT} = flow_in(xylem.auxil);

flow_in(flow::XylemHydraulicsAuxilNSS{FT}) where {FT} = flow.flow[1];

flow_in(flow::XylemHydraulicsAuxilSS{FT}) where {FT} = flow.flow;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-22: add function to read the flow rate out of the xylem
#
#######################################################################################################################################################################################################
"""

    flow_out(organ)

Return the flow rate out of the organ

"""
function flow_out end;

flow_out(xylem::XylemHydraulics{FT}) where {FT} = flow_out(xylem.auxil);

flow_out(x_aux::XylemHydraulicsAuxilNSS{FT}) where {FT} = x_aux.flow[end];

flow_out(x_aux::XylemHydraulicsAuxilSS{FT}) where {FT} = x_aux.flow;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function to set the flow rate profile in the xylem based on the flow rate out of the xylem
#     2024-Feb-28: add LAI <= 0 control
#
#######################################################################################################################################################################################################
"""

    set_flow_profile!(xylem::XylemHydraulics{FT}, flow::FT) where {FT}

Set the flow rate profile in the xylem, given
- `xylem` `XylemHydraulics` type struct
- `flow` Flow rate out of the xylem

"""
function set_flow_profile! end;

set_flow_profile!(xylem::XylemHydraulics{FT}, flow::FT) where {FT} = (
    if xylem.trait.area <= 0
        return nothing
    end;

    # update the flow profile calculation only if xylem area > 0
    set_flow_profile!(xylem.auxil, flow);

    return nothing
);

set_flow_profile!(x_aux::XylemHydraulicsAuxilNSS{FT}, flow::FT) where {FT} = (
    x_aux.flow[end] = flow;
    for i in length(x_aux.flow_buffer):-1:1
        x_aux.flow[i] = x_aux.flow[i+1] - x_aux.flow_buffer[i];
    end;

    return nothing
);

set_flow_profile!(x_aux::XylemHydraulicsAuxilSS{FT}, flow::FT) where {FT} = (
    x_aux.flow = flow;

    return nothing
);

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: add function from xylem_flow_profile! to set_leaf_flow! to be more specific in function name
#
#######################################################################################################################################################################################################
"""

    set_flow!(mode::SteadyStateFlow, f_out::FT) where {FT<:AbstractFloat}
    set_flow!(mode::NonSteadyStateFlow, f_out::FT) where {FT<:AbstractFloat}

Set flow profile based on the flow mode, given
- `mode` `SteadyStateFlow` or `NonSteadyStateFlow` type struct
- `f_out` flow rate exiting the organ

"""
function set_flow_out! end

set_flow_out!(mode::SteadyStateFlow, f_out::FT) where {FT<:AbstractFloat} = (mode.flow = f_out; return nothing);

set_flow_out!(mode::NonSteadyStateFlow, f_out::FT) where {FT<:AbstractFloat} = (mode.f_out = f_out; return nothing);

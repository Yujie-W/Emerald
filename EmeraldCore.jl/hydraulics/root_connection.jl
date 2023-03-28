#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Mar-28: add function to disconnect root from soil
#
#######################################################################################################################################################################################################
"""

    disconnect!(root::Root{FT}) where {FT<:AbstractFloat}

Disconnect root from soil, given
- `root` `Root` type struct

"""
function disconnect! end

disconnect!(root::Root{FT}) where {FT<:AbstractFloat} = (
    root._isconnected = false;
    disconnect!(root.HS, root.HS.FLOW);

    return nothing
);

disconnect!(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = (
    # update the pressure
    hs._p_element .= hs._p_storage;

    # update the flow
    mode.f_in = 0;
    mode.f_out = 0;
    mode._f_buffer .= 0;
    mode._f_element .= 0;
    mode._f_sum .= 0;

    return nothing
);

disconnect!(hs::RootHydraulics{FT}, mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = (
    mode.flow = 0;

    return nothing
);

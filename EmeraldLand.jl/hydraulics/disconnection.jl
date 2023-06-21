#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Mar-28: add function to disconnect root from soil
#     2023-Mar-28: add function to disconnect others (by setting flows to 0)
#
#######################################################################################################################################################################################################
"""

    disconnect!(organ::Union{Leaf{FT},Leaves2D{FT},Stem{FT}}) where {FT<:AbstractFloat}
    disconnect!(organ::Leaves1D{FT}) where {FT<:AbstractFloat}
    disconnect!(organ::Root{FT}) where {FT<:AbstractFloat}

Disconnect root from soil (and set othes' flow to 0), given
- `organ` Root, stem, or leaf

"""
function disconnect! end

disconnect!(organ::Union{Leaf{FT},Leaves2D{FT},Stem{FT}}) where {FT<:AbstractFloat} = (
    disconnect!(organ.HS, organ.HS.FLOW);

    return nothing
);

disconnect!(organ::Leaves1D{FT}) where {FT<:AbstractFloat} = (
    disconnect!(organ.HS1, organ.HS1.FLOW);
    disconnect!(organ.HS2, organ.HS2.FLOW);

    return nothing
);

disconnect!(organ::Root{FT}) where {FT<:AbstractFloat} = (
    organ._isconnected = false;
    disconnect!(organ.HS, organ.HS.FLOW);

    return nothing
);

disconnect!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = (
    # update the pressure
    hs.p_leaf = hs._p_storage;

    # update the flow
    mode.f_in = 0;
    mode.f_out = 0;
    mode._f_buffer .= 0;
    mode._f_element .= 0;
    mode._f_sum .= 0;

    return nothing
);

disconnect!(hs::Union{RootHydraulics{FT}, StemHydraulics{FT}}, mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = (
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

disconnect!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = (
    mode.flow = 0;

    return nothing
);

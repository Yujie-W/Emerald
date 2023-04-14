#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Mar-28: add function to disconnect root from soil
#     2023-Mar-28: add function to disconnect others (by setting flows to 0)
#
#######################################################################################################################################################################################################
"""

    disconnect!(organ::Union{Leaf{FT},Leaves2D{FT},Stem{FT}}) where {FT}
    disconnect!(organ::Leaves1D{FT}) where {FT}
    disconnect!(organ::Root{FT}) where {FT}

Disconnect root from soil (and set othes' flow to 0), given
- `organ` Root, stem, or leaf

"""
function disconnect! end

disconnect!(organ::Union{Leaf{FT},Leaves2D{FT},Stem{FT}}) where {FT} = (
    disconnect!(organ.HS, organ.HS.FLOW);

    return nothing
);

disconnect!(organ::Leaves1D{FT}) where {FT} = (
    disconnect!(organ.HS1, organ.HS1.FLOW);
    disconnect!(organ.HS2, organ.HS2.FLOW);

    return nothing
);

disconnect!(organ::Root{FT}) where {FT} = (
    organ._isconnected = false;
    disconnect!(organ.HS, organ.HS.FLOW);

    return nothing
);

disconnect!(hs::LeafHydraulics{FT,DIM_XYLEM}, mode::NonSteadyStateFlow{FT,1}) where {FT,DIM_XYLEM} = (
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

disconnect!(hs::Union{RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}, mode::NonSteadyStateFlow{FT,DIM_XYLEM}) where {FT,DIM_XYLEM} = (
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

disconnect!(hs::Union{LeafHydraulics{FT,DIM_XYLEM}, RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}, mode::SteadyStateFlow{FT}) where {FT,DIM_XYLEM} = (
    mode.flow = 0;

    return nothing
);

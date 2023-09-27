#=

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: add function to extract flow rate
#     2022-May-31: add method to extract flow rate from steady state flow to use with upstream flow
#     2022-May-31: add method to extract flow rate from non-steady state flow to use with upstream flow
#     2022-May-31: add method to extract flow rate from hydraulic system
#     2022-May-31: add method to extract flow rate from organ
#     2022-May-31: add method to extract flow rate from leaves
#     2022-May-31: add method to extract flow rate from branches
#     2022-Jun-30: add method to extract flow rate from Leaves1D
#     2022-Jun-30: add support for Leaves2D
#     2022-Jun-30: rename Leaf to Leaves2D to support ML*SPAC
#     2022-Jul-08: deflate documentations
#     2022-Jul-15: rename xylem_flow to flow_in to be more descriptive
#
#######################################################################################################################################################################################################
"""

    flow_in(organ::Union{Leaf2{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT}
    flow_in(organ::Leaves1D{FT}) where {FT}
    flow_in(organs::Vector{Leaves2D{FT}}) where {FT}
    flow_in(organs::Vector{Stem{FT}}) where {FT}

Return the flow rate, given
- `organ` `Leaf2`, `Leaves1D`, `Leaves2D`, `Root`, or `Stem` type struct
- `organs` Vector of `Leaves2D` or `Stem` type struct

"""
function flow_in end

flow_in(organ::Union{Leaf2{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT} = flow_in(organ.HS);

flow_in(organ::Leaves1D{FT}) where {FT} = (flow_in(organ.HS), flow_in(organ.HS2));

flow_in(organs::Vector{Leaves2D{FT}}) where {FT} = (
    f_sum::FT = 0;
    for i in eachindex(organs)
        f_sum += flow_in(organs[i]) * organs[i].HS.AREA;
    end;

    return f_sum
);

flow_in(organs::Vector{Stem{FT}}) where {FT} = (
    f_sum::FT = 0;
    for i in eachindex(organs)
        f_sum += flow_in(organs[i]);
    end;

    return f_sum
);

flow_in(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}) where {FT} = flow_in(hs.FLOW);

flow_in(mode::SteadyStateFlow{FT}) where {FT} = mode.flow;

flow_in(mode::NonSteadyStateFlow{FT}) where {FT} = mode.f_in;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-15: add function to read water exiting the leaf
#     2022-Jul-15: rename to flow_out to be more descriptive
#
#######################################################################################################################################################################################################
"""

    flow_out(lf::Union{Leaf2{FT}, Leaves2D{FT}}) where {FT}

Return the net flow that escape from the leaf, given
- `lf` `Leaf2`, `Leaves2D`, `Root`, or `Stem` type organ

"""
function flow_out end

flow_out(organ::Union{Leaf2{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT} = flow_out(organ.HS.FLOW);

flow_out(mode::SteadyStateFlow{FT}) where {FT} = mode.flow;

flow_out(mode::NonSteadyStateFlow{FT}) where {FT} = mode.f_out;

=#

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: rename method to leaf_flow_profile! to be more specific in function name
#
#######################################################################################################################################################################################################
"""

    set_stem_flow_out!(spac::MonoElementSPAC{FT}) where {FT}
    set_stem_flow_out!(spac::MultiLayerSPAC{FT}) where {FT}

Set the flow out from each stem, given
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type struct

"""
function set_stem_flow_out! end

set_stem_flow_out!(spac::MonoElementSPAC{FT}) where {FT} = (
    (; LEAF, STEM) = spac;

    set_flow_out!(STEM.HS.FLOW, flow_in(LEAF) * LEAF.HS.AREA);

    return nothing
);

set_stem_flow_out!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; BRANCHES, LEAVES) = spac;

    for _i in eachindex(LEAVES)
        set_flow_out!(BRANCHES[_i].HS.FLOW, flow_in(LEAVES[_i]) * LEAVES[_i].HS.AREA);
    end;

    return nothing
);

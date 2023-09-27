#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: rename method to leaf_flow_profile! to be more specific in function name
#     2023-Sep-12: set stem flow out to 0 if there is no leaves
#
#######################################################################################################################################################################################################
"""

    set_stem_flow_out!(spac::MultiLayerSPAC{FT}) where {FT}

Set the flow out from each stem, given
- `spac` `MultiLayerSPAC` type struct

"""
function set_stem_flow_out! end

set_stem_flow_out!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; BRANCHES, CANOPY, LEAVES) = spac;

    if CANOPY.lai > 0
        for _i in eachindex(BRANCHES)
            set_flow_out!(BRANCHES[_i].HS.FLOW, flow_in(LEAVES[_i]) * LEAVES[_i].HS.AREA);
        end;
    else
        for _i in eachindex(BRANCHES)
            set_flow_out!(BRANCHES[_i].HS.FLOW, FT(0));
        end;
    end;

    return nothing
);

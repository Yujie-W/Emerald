# This file contains functions related to stem flow profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-05: add function stem_flow_profiles!
#     2023-Sep-28: update the branch and trunk flow profile at the same time
#
#######################################################################################################################################################################################################
"""

    stem_flow_profiles!(spac::MultiLayerSPAC{FT}) where {FT}

Set up stem flow profile, given
- `spac` `MultiLayerSPAC` type struct

"""
function stem_flow_profiles!(spac::MultiLayerSPAC{FT}) where {FT}
    (; BRANCHES, LEAVES, TRUNK) = spac;

    sum_f::FT = 0;
    for i in eachindex(BRANCHES)
        set_flow_profile!(BRANCHES[i].NS.xylem, flow_in(LEAVES[i].NS.xylem));
        sum_f += flow_out(BRANCHES[i].NS.xylem);
    end;

    set_flow_profile!(TRUNK.NS.xylem, sum_f);

    return nothing
end;

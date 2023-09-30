# This file contains functions related to stem flow profile

flow_in(stem::Stem{FT}) where {FT} = flow_in(stem.xylem);

flow_in(stem::Stem2{FT}) where {FT} = flow_in(stem.NS);

flow_out(stem::Stem{FT}) where {FT} = flow_out(stem.xylem);

flow_out(stem::Stem2{FT}) where {FT} = flow_out(stem.NS);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-05: add function stem_flow_profiles!
#     2023-Sep-28: update the branch and trunk flow profile at the same time
#     2023-Sep-30: subtract the flow into TRUNK from the JUNCTION.auxil.∂w∂t
#
#######################################################################################################################################################################################################
"""

    stem_flow_profiles!(spac::MultiLayerSPAC{FT}) where {FT}

Set up stem flow profile, given
- `spac` `MultiLayerSPAC` type struct

"""
function stem_flow_profiles!(spac::MultiLayerSPAC{FT}) where {FT}
    (; BRANCHES, JUNCTION, LEAVES, TRUNK) = spac;

    sum_f::FT = 0;
    for i in eachindex(BRANCHES)
        set_flow_profile!(BRANCHES[i].NS.xylem, flow_in(LEAVES[i]));
        sum_f += flow_out(BRANCHES[i]);
    end;

    set_flow_profile!(TRUNK.NS.xylem, sum_f);
    JUNCTION.auxil.∂w∂t -= flow_in(TRUNK);

    return nothing
end;

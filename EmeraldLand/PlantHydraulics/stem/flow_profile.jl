# This file contains functions related to stem flow profile

flow_in(stem::Stem{FT}) where {FT} = flow_in(stem.xylem);

flow_out(stem::Stem{FT}) where {FT} = flow_out(stem.xylem);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-05: add function stem_flow_profiles!
#     2023-Sep-28: update the branch and trunk flow profile at the same time
#     2023-Sep-30: subtract the flow into trunk from the junction.auxil.∂w∂t
#     2024-Feb-28: add LAI <= 0 control
#
#######################################################################################################################################################################################################
"""

    stem_flow_profiles!(spac::BulkSPAC{FT}) where {FT}

Set up stem flow profile, given
- `spac` `BulkSPAC` type struct

"""
function stem_flow_profiles!(spac::BulkSPAC{FT}) where {FT}
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    # run the flow profile calculation for each stem layer only if LAI > 0
    branches = spac.plant.branches;
    junction = spac.plant.junction;
    leaves = spac.plant.leaves;
    trunk = spac.plant.trunk;

    sum_f::FT = 0;
    for i in eachindex(branches)
        set_flow_profile!(branches[i].xylem, flow_in(leaves[i]));
        sum_f += flow_out(branches[i]);
    end;

    set_flow_profile!((trunk).xylem, sum_f);
    junction.auxil.∂w∂t -= flow_in(trunk);

    return nothing
end;

# This file contains functions to update the water budget of the root-trunk junction

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function junction_water_budget!
#
#######################################################################################################################################################################################################
"""

    junction_water_budget!(spac::BulkSPAC{FT}) where {FT}

Update the water budget of the root-trunk junction, given
- `spac` `BulkSPAC` type struct

"""
function junction_water_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    # compute the total flow into the junction as the sum of root flow out minus the trunk flow in
    sum_q::FT = 0;
    for root in spac.plant.roots
        sum_q += flow_out(root) * δt;
    end;
    sum_q -= flow_in(spac.plant.trunk) * δt;

    spac.plant.junction.state.v_storage += sum_q;

    return nothing
end;

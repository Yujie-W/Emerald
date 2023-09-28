# This file contains functions to update the water budget of the root-trunk junction

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function junction_water_budget!
#
#######################################################################################################################################################################################################
"""

    junction_water_budget!(spac::MultiLayerSPAC{FT}) where {FT}

Update the water budget of the root-trunk junction, given
- `spac` `MultiLayerSPAC` type struct

"""
function junction_water_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    (; JUNCTION, ROOTS, TRUNK) = spac;

    # compute the total flow into the junction as the sum of root flow out minus the trunk flow in
    sum_q::FT = 0;
    for root in ROOTS
        sum_q += flow_out(root.NS.xylem) * δt;
    end;
    sum_q -= flow_in(TRUNK.NS.xylem) * δt;

    JUNCTION.state.v_storage += sum_q;
    JUNCTION.auxil.pressure = capacitance_pressure(JUNCTION.state.pv, JUNCTION.state.v_storage / JUNCTION.state.v_max, JUNCTION.auxil.t);

    return nothing
end;

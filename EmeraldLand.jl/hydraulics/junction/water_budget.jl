# This file contains functions to update the water budget of the root-trunk junction


function junction_water_budget!(spac::MultiLayerSPAC{FT}) where {FT}
    (; JUNCTION, ROOTS, TRUNK) = spac;

    # compute the total flow into the junction as the sum of root flow out minus the trunk flow in
    sum_q::FT = 0;
    for root in ROOTS
        sum_q += root.∫∂w∂t_out;
    end;
    sum_q -= TRUNK.∫∂w∂t_in;

    JUNCTION.state.v_storage += sum_q;
    JUNCTION.auxil.pressure = capacitance_pressure(JUNCTION.state.pv, JUNCTION.state.v_storage / JUNCTION.state.v_max, JUNCTION.auxil.t);

    return nothing
end;

# THis file contains functions to update the water budget of the leaf

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add leaf_water_budget! function
#
#######################################################################################################################################################################################################
"""

    leaf_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}

Set the flow profile of the leaf, given
- `spac` `MultiLayerSPAC` type struct
- `Δt` time step

"""
function leaf_water_budget! end

leaf_water_budget!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT} = (
    if spac.CANOPY.lai > 0
        # do this way to avoid memory allocation of a [nothing...] vector
        for leaf in spac.LEAVES
            leaf_water_budget!(leaf, (leaf).xylem.auxil, Δt);
        end;
    end;

    return nothing
);

leaf_water_budget!(leaf::Leaf{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, Δt::FT) where {FT} = (
    f_vis = relative_viscosity(leaf.energy.auxil.t);

    # make sure the buffer rate does not drain or overflow the capacictance
    # TODO: add this to time_stepper! function, otherwise the water budget will not be consvered
    if leaf.capacitor.auxil.flow > 0 && leaf.capacitor.state.v_storage * leaf.xylem.state.area <= leaf.capacitor.auxil.flow * Δt
        @warn "The capacitance buffer is drained, use only half of the remaining water in the buffer!";
        leaf.capacitor.auxil.flow = (leaf.capacitor.state.v_storage * leaf.xylem.state.area / 2) / Δt;
    end;

    # update the integrators of the flow
    leaf.∫∂w∂t_in += flow_in(x_aux) * Δt;
    leaf.∫∂w∂t_out += (flow_out(x_aux) + leaf.capacitor.auxil.flow) * Δt;

    # update storage and the tissue pressure (p_storage)
    leaf.capacitor.state.v_storage -= leaf.capacitor.auxil.flow * Δt / leaf.xylem.state.area;
    leaf.capacitor.state.p = capacitance_pressure(leaf.capacitor.state.pv, leaf.capacitor.state.v_storage / leaf.capacitor.state.v_max, leaf.energy.auxil.t);
    leaf.capacitor.auxil.flow = (x_aux.pressure[end] - leaf.capacitor.state.p) * leaf.capacitor.state.pv.k_refill / f_vis * leaf.capacitor.state.v_max * leaf.xylem.state.area;

    return nothing
);

leaf_water_budget!(leaf::Leaf{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, Δt::FT) where {FT} = (
    leaf.∫∂w∂t_in += flow_in(x_aux) * Δt;
    leaf.∫∂w∂t_out += flow_out(x_aux) * Δt;

    return nothing
);

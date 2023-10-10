# THis file contains functions to update the water budget of the leaf

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-27: add leaf_water_budget! function
#     2023-Sep-28: make the buffer flow to be flow out from the capacitor
#
#######################################################################################################################################################################################################
"""

    leaf_water_budget!(leaf::Leaf{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, δt::FT) where {FT}
    leaf_water_budget!(leaf::Leaf{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, δt::FT) where {FT}

Set the flow profile of the leaf, given
- `leaf` `Leaf` type struct
- `x_aux` `XylemHydraulicsAuxilNSS` or `XylemHydraulicsAuxilSS` type struct
- `δt` time step

"""
function leaf_water_budget! end;

leaf_water_budget!(leaf::Leaf{FT}, x_aux::XylemHydraulicsAuxilNSS{FT}, δt::FT) where {FT} = (
    # make sure the buffer rate does not drain or overflow the capacictance
    # TODO: add this to time_stepper! function, otherwise the water budget will not be consvered
    if leaf.capacitor.auxil.flow > 0 && leaf.capacitor.state.v_storage * leaf.xylem.state.area <= leaf.capacitor.auxil.flow * δt
        @warn "The capacitance buffer is drained, use only half of the remaining water in the buffer!";
        leaf.capacitor.auxil.flow = (leaf.capacitor.state.v_storage * leaf.xylem.state.area / 2) / δt;
    end;

    # update the integrators of the flow
    leaf.flux.auxil.∫∂w∂t_out += flow_out(leaf) * δt;

    # update storage and the tissue pressure (p_storage)
    leaf.capacitor.state.v_storage -= leaf.capacitor.auxil.flow * δt / leaf.xylem.state.area;

    return nothing
);

leaf_water_budget!(leaf::Leaf{FT}, x_aux::XylemHydraulicsAuxilSS{FT}, δt::FT) where {FT} = (
    leaf.flux.auxil.∫∂w∂t_out += flow_out(leaf) * δt;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add leaf_water_budgets! function
#
#######################################################################################################################################################################################################
"""

    leaf_water_budgets!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Set the flow profile of each leaf, given
- `spac` `MultiLayerSPAC` type struct
- `δt` time step

"""
function leaf_water_budgets!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    if spac.CANOPY.lai > 0
        # do this way to avoid memory allocation of a [nothing...] vector
        for leaf in spac.LEAVES
            leaf_water_budget!(leaf, (leaf).xylem.auxil, δt);
        end;
    end;

    return nothing
end;

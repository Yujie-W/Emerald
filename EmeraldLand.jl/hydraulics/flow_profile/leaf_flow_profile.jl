#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: rename method to leaf_flow_profile! to be more specific in function name
#     2023-Sep-12: update leaf flow profile only if lai > 0
#
#######################################################################################################################################################################################################
"""

    leaf_flow_profile!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
    leaf_flow_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}

Set the flow profile of the leaf, given
- `config` `SPACConfiguration` type struct
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type struct
- `Δt` time step

"""
function leaf_flow_profile! end

leaf_flow_profile!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}, Δt::FT) where {FT} = (
    set_leaf_flow_out!(config, spac);
    leaf_flow_profile!(spac.LEAF, Δt);

    return nothing
);

leaf_flow_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, Δt::FT) where {FT} = (
    if spac.CANOPY.lai > 0
        set_leaf_flow_out!(config, spac);
        leaf_flow_profile!.(spac.LEAVES, Δt);
    end;

    return nothing
);

leaf_flow_profile!(organ::Leaf{FT}, Δt::FT) where {FT<:AbstractFloat} = leaf_flow_profile!(organ.HS, organ.t, Δt);

leaf_flow_profile!(organ::Leaves1D{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    (; HS, HS2) = organ;

    leaf_flow_profile!(HS, organ.t[1], Δt);
    leaf_flow_profile!(HS2, organ.t[2], Δt);

    return nothing
);

leaf_flow_profile!(organ::Leaves2D{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    leaf_flow_profile!(organ.HS, organ.t, Δt);
    organ.∫∂w∂t_in += flow_in(organ) * Δt;
    organ.∫∂w∂t_out += flow_out(organ) * Δt;

    return nothing
);

leaf_flow_profile!(hs::LeafHydraulics{FT}, t::FT, Δt::FT) where {FT<:AbstractFloat} = leaf_flow_profile!(hs, hs.FLOW, t, Δt);

leaf_flow_profile!(hs::LeafHydraulics{FT}, mode::SteadyStateFlow{FT}, t::FT, Δt::FT) where {FT<:AbstractFloat} = nothing;

leaf_flow_profile!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}, t::FT, Δt::FT) where {FT<:AbstractFloat} = (
    (; PVC, V_MAXIMUM) = hs;

    _f_vis = relative_viscosity(t);

    # compute the flow rate from capacitance buffer
    mode._f_buffer[1] = (hs._p_storage - hs.p_leaf) * capacitance_buffer(PVC) / _f_vis * V_MAXIMUM;

    # make sure the buffer rate does not drain or overflow the capacictance
    if (mode._f_buffer[1] > 0) && (hs.v_storage <= mode._f_buffer[1] * Δt)
        mode._f_buffer[1] = (hs.v_storage - eps(FT)) / Δt;
    end;

    # update storage and the tissue pressure (p_storage)
    hs.v_storage -= mode._f_buffer[1] * Δt;
    hs._p_storage = xylem_pressure(PVC, hs.v_storage/V_MAXIMUM, t);

    # update flow into the tissue an along the xylem (same as f_in)
    mode.f_in = mode.f_out - mode._f_buffer[1];
    mode._f_element .= mode.f_in;

    return nothing
);

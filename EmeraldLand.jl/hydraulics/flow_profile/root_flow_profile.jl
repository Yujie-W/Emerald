#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: rename method to leaf_flow_profile! to be more specific in function name
#
#######################################################################################################################################################################################################
"""

    set_root_flow_out!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}) where {FT}
    set_root_flow_out!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Set the flow profile of each root, given
- `config` `SPACConfiguration` type struct
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type struct

"""
function root_flow_profile! end

root_flow_profile!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}, Δt::FT) where {FT} = (
    (; ROOT) = spac;

    set_root_flow_out!(spac);
    root_flow_profile!(ROOT, Δt);

    return nothing
);

root_flow_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, Δt::FT) where {FT} = (
    (; ROOTS) = spac;

    set_root_flow_out!(config, spac);
    root_flow_profile!.(ROOTS, Δt);

    return nothing
);

root_flow_profile!(organ::Root{FT}, Δt::FT) where {FT} = (
    if organ._isconnected
        root_flow_profile!(organ.HS, organ.t, Δt);
        organ.∫∂w∂t_in += flow_in(organ) * Δt;
        organ.∫∂w∂t_out += flow_out(organ) * Δt;
    end;

    return nothing
);

root_flow_profile!(hs::RootHydraulics{FT}, T::FT, Δt::FT) where {FT} = root_flow_profile!(hs, hs.FLOW, T, Δt);

root_flow_profile!(hs::RootHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT, Δt::FT) where {FT} = nothing;

root_flow_profile!(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT, Δt::FT) where {FT} = (
    (; DIM_XYLEM, PVC, V_MAXIMUM) = hs;

    _f_vis = relative_viscosity(T);

    # update storage volume and pressure per slice
    _f_sum::FT = 0;
    for _i in DIM_XYLEM:-1:1
        mode._f_buffer[_i] = (hs._p_storage[_i] - hs._p_element[_i]) * capacitance_buffer(PVC) / _f_vis * V_MAXIMUM[_i];

        # make sure the buffer rate does not drain or overflow the capacictance
        if (mode._f_buffer[_i] > 0) && (hs.v_storage[_i] <= mode._f_buffer[_i] * Δt)
            mode._f_buffer[_i] = hs.v_storage[_i] / Δt;
        end;

        mode._f_sum[_i] = _f_sum;
        hs.v_storage[_i] -= mode._f_buffer[_i] * Δt;
        hs._p_storage[_i] = xylem_pressure(PVC, hs.v_storage[_i]/V_MAXIMUM[_i], T);
        mode._f_element[_i] = mode.f_out - _f_sum;
        _f_sum += mode._f_buffer[_i];
    end;

    # update flow into the tissue
    mode.f_in = mode.f_out - _f_sum;

    return nothing
);



function stem_flow_profile! end

stem_flow_profile!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT} = (
    set_stem_flow_out!(spac);
    stem_flow_profile!(spac.STEM, Δt);

    return nothing
);

stem_flow_profile!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT} = (
    (; BRANCHES, TRUNK) = spac;

    set_stem_flow_out!(spac);
    stem_flow_profile!.(BRANCHES, Δt);

    # update the flow profile in the trunk as well
    set_flow_out!(TRUNK.HS.FLOW, flow_in(BRANCHES));
    stem_flow_profile!(TRUNK, Δt);

    return nothing
);

stem_flow_profile!(organ::Stem{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    stem_flow_profile!(organ.HS, organ.t, Δt);
    organ.∫∂w∂t_in += flow_in(organ) * Δt;
    organ.∫∂w∂t_out += flow_out(organ) * Δt;

    return nothing
);

stem_flow_profile!(hs::StemHydraulics{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat} = stem_flow_profile!(hs, hs.FLOW, T, Δt);

stem_flow_profile!(hs::StemHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat} = nothing;

stem_flow_profile!(hs::StemHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat} = (
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

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: rename method to leaf_flow_profile! to be more specific in function name
#     2023-Sep-12: redo the root water function by disconnecting the roots with lowest conductance during a drought
#
#######################################################################################################################################################################################################
"""

    set_root_flow_out!(spac::MonoElementSPAC{FT}) where {FT}
    set_root_flow_out!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Set the flow out from each root, given
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type struct
- `config` `SPACConfiguration` type struct

"""
function set_root_flow_out! end

set_root_flow_out!(spac::MonoElementSPAC{FT}) where {FT} = (
    (; STEM, ROOT) = spac;

    set_flow_out!(ROOT.HS.FLOW, flow_in(STEM));

    return nothing
);

set_root_flow_out!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    # very first step here: if soil is too dry, disconnect root from soil
    disconnect_roots!(config, spac);
    disconnect_spac!(spac);
    if !spac._root_connection
        return nothing
    end;

    # if spac is still connected to soil, continue
    (; ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;
    _f_sum = flow_in(TRUNK);

    # update root buffer rates to get an initial guess (flow rate not changing now as time step is set to 0)
    root_flow_profile!.(ROOTS, FT(0));

    # recalculate the flow profiles to make sure sum are the same as f_sum
    # _use_second_solver = false;
    _count = 0;
    while _count < 20
        _count += 1;
        # sync the values to ks, ps, and qs
        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
            if _root._isconnected
                set_flow_out!(_root.HS.FLOW, spac._fs[_i]);
                if _root.HS.FLOW isa NonSteadyStateFlow
                    _root.HS.FLOW._f_element .= _root.HS.FLOW.f_out .- _root.HS.FLOW._f_sum;
                end;
            else
                spac._fs[_i] = 0;
            end;
            spac._ps[_i],spac._ks[_i] = root_pk(_root, _slayer);
        end;

        # use ps and ks to compute the Δf to adjust
        _pm = nanmean(spac._ps);
        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            if _root._isconnected
                spac._fs[_i] -= (_pm - spac._ps[_i]) * spac._ks[_i];
            else
                spac._fs[_i] = 0;
            end;
        end;

        # adjust the fs so that sum(fs) = f_sum
        _f_diff = sum(spac._fs) - _f_sum;
        if (abs(_f_diff) < FT(1e-6)) && (nanmax(spac._ps) - nanmin(spac._ps) < 1e-4)
            break
        end;
        _k_sum  = sum(spac._ks);
        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            if _root._isconnected
                spac._fs[_i] -= _f_diff * spac._ks[1] / _k_sum;
            else
                spac._fs[_i] = 0;
            end;
        end;

        # if the flow and pressure profiles cannot converge, disconnect the soil with lowest conductance
        if _count > 15
            _root_to_disconnect = 0;
            _min_kr = 1;
            for _i in eachindex(ROOTS)
                _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
                _kr = relative_hydraulic_conductance(_slayer.VC, _slayer.θ);
                if _kr < _min_kr
                    _min_kr = _kr;
                    _root_to_disconnect = _i;
                end;
            end;
            # disconnect the root with lowest conductance, determine if the spac is still connected to soil
            disconnect!(ROOTS[_root_to_disconnect]);
            disconnect_spac!(spac);

            # if no roots are connected to soil, stop the while loop
            if spac._root_connection
                break
            end;

            # restart the counting
            _count = 0;
        end;

        # if _count == 20
        #     @warn "Root flow profile solver did not converge, !";
        #     _use_second_solver = true;
        # end;
    end;

    #=
    # disable the second solver for now
    # use second solver to solve for the flow rates (when SWC differs a lot among layers)
    if _use_second_solver
        @inline diff_p_root(ind::Int, e::FT, p::FT) where {FT<:AbstractFloat} = (
            _root = ROOTS[ind];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[ind]];
            if _root._isconnected
                set_flow_out!(_root.HS.FLOW, e);
                _root.HS.FLOW._f_element .= _root.HS.FLOW.f_out .- _root.HS.FLOW._f_sum;
            end;
            (_p,_) = root_pk(_root, _slayer);

            return _p - p
        );

        @inline diff_e_root(p::FT) where {FT<:AbstractFloat} = (
            _sum::FT = 0;
            for _i in eachindex(ROOTS_INDEX)
                _root = ROOTS[_i];
                if _root._isconnected
                    _f(e) = diff_p_root(_i, e, p);
                    _tol = SolutionTolerance{FT}(1e-8, 50);
                    _met = NewtonBisectionMethod{FT}(x_min = -1000, x_max = 1000, x_ini = 0);
                    _sol = find_zero(_f, _met, _tol);
                    _sum += _sol;
                end;
            end;

            return _sum - _f_sum
        );

        _tol = SolutionTolerance{FT}(1e-8, 50);
        _met = NewtonBisectionMethod{FT}(x_min = -1000, x_max = 1000, x_ini = 0);
        _p_r = find_zero(diff_e_root, _met, _tol);

        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];

            if _root._isconnected
                _f(e) = diff_p_root(_i, e, _p_r);
                _tol = SolutionTolerance{FT}(1e-8, 50);
                _met = NewtonBisectionMethod{FT}(x_min = -1000, x_max = 1000, x_ini = 0);
                spac._fs[_i] = find_zero(_f, _met, _tol);
                set_flow_out!(_root.HS.FLOW, spac._fs[_i]);
                _root.HS.FLOW._f_element .= _root.HS.FLOW.f_out .- _root.HS.FLOW._f_sum;
            else
                spac._fs[_i] = 0;
            end;

            spac._ps[_i],spac._ks[_i] = root_pk(_root, _slayer);
        end;
    end;
    =#

    # update root buffer rates again
    if spac._root_connection
        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            if _root._isconnected
                set_flow_out!(ROOTS[_i].HS.FLOW, spac._fs[_i]);
            end;
        end;
    end;

    return nothing
);

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-11: rename method to leaf_flow_profile! to be more specific in function name
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
    (; KR_THRESHOLD) = config;
    (; BRANCHES, LEAVES, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    _f_sum = flow_in(TRUNK);

    # very first step here: if soil is too dry, disconnect root from soil
    _connected = 0;
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        _ψ_soil = soil_ψ_25(_slayer.VC, _slayer.θ) * relative_surface_tension(_slayer.t);
        _p_crit = xylem_pressure(_root.HS.VC, KR_THRESHOLD) * relative_surface_tension(_root.t);
        if _ψ_soil <= _p_crit
            disconnect!(_root);
        else
            _root._isconnected = true;
            _connected += 1;
        end;
    end;

    # if all roots are disconnected, set all flows to 0
    # TODO: if leaf shedding is allowed, remember to update leaf area when roots are reconnected to the soil
    if _connected > 0
        if !spac._root_connection
            @warn "Roots are now reconnected to soil!";
        end;
        spac._root_connection = true;
    else
        @warn "Roots are now all disconnected from soil!"
        spac._root_connection = false;
        disconnect!(TRUNK);
        disconnect!.(BRANCHES);
        disconnect!.(LEAVES);

        return nothing
    end;

    # update root buffer rates to get an initial guess (flow rate not changing now as time step is set to 0)
    root_flow_profile!.(ROOTS, FT(0));

    # recalculate the flow profiles to make sure sum are the same as f_sum
    _use_second_solver = false;
    for _count in 1:20
        # sync the values to ks, ps, and qs
        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
            if _root._isconnected
                set_flow_out!(_root.HS.FLOW, spac._fs[_i]);
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

        if _count == 20
            _use_second_solver = true
        end;
    end;

    # use second solver to solve for the flow rates (when SWC differs a lot among layers)
    if _use_second_solver
        @inline diff_p_root(ind::Int, e::FT, p::FT) where {FT<:AbstractFloat} = (
            _root = ROOTS[ind];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[ind]];
            if _root._isconnected
                set_flow_out!(ROOTS[ind].HS.FLOW, e);
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
                set_flow_out!(ROOTS[_i].HS.FLOW, spac._fs[_i]);
            else
                spac._fs[_i] = 0;
            end;

            spac._ps[_i],spac._ks[_i] = root_pk(_root, _slayer);
        end;
    end;

    # update root buffer rates again
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        if _root._isconnected
            set_flow_out!(ROOTS[_i].HS.FLOW, spac._fs[_i]);
        end;
    end;

    return nothing
);

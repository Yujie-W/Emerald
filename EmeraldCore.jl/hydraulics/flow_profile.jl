#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: add function to extract flow rate
#     2022-May-31: add method to extract flow rate from steady state flow to use with upstream flow
#     2022-May-31: add method to extract flow rate from non-steady state flow to use with upstream flow
#     2022-May-31: add method to extract flow rate from hydraulic system
#     2022-May-31: add method to extract flow rate from organ
#     2022-May-31: add method to extract flow rate from leaves
#     2022-May-31: add method to extract flow rate from branches
#     2022-Jun-30: add method to extract flow rate from Leaves1D
#     2022-Jun-30: add support for Leaves2D
#     2022-Jun-30: rename Leaf to Leaves2D to support ML*SPAC
#     2022-Jul-08: deflate documentations
#     2022-Jul-15: rename xylem_flow to flow_in to be more descriptive
#
#######################################################################################################################################################################################################
"""

    flow_in(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT}
    flow_in(organ::Leaves1D{FT}) where {FT}
    flow_in(organs::Vector{Leaves2D{FT}}) where {FT}
    flow_in(organs::Vector{Stem{FT}}) where {FT}

Return the flow rate, given
- `organ` `Leaf`, `Leaves1D`, `Leaves2D`, `Root`, or `Stem` type struct
- `organs` Vector of `Leaves2D` or `Stem` type struct

"""
function flow_in end

flow_in(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT} = flow_in(organ.HS);

flow_in(organ::Leaves1D{FT}) where {FT} = (flow_in(organ.HS), flow_in(organ.HS2));

flow_in(organs::Vector{Leaves2D{FT}}) where {FT} = (
    _f_sum::FT = 0;
    for _i in eachindex(organs)
        _f_sum += flow_in(organs[_i]) * organs[_i].HS.AREA;
    end;

    return _f_sum
);

flow_in(organs::Vector{Stem{FT}}) where {FT} = (
    _f_sum::FT = 0;
    for _i in eachindex(organs)
        _f_sum += flow_in(organs[_i]);
    end;

    return _f_sum
);

flow_in(hs::Union{LeafHydraulics{FT,DIM_XYLEM}, RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}) where {FT,DIM_XYLEM} = flow_in(hs.FLOW);

flow_in(mode::SteadyStateFlow{FT}) where {FT} = mode.flow;

flow_in(mode::NonSteadyStateFlow{FT,DIM_XYLEM}) where {FT,DIM_XYLEM} = mode.f_in;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-15: add function to read water exiting the leaf
#     2022-Jul-15: rename to flow_out to be more descriptive
#
#######################################################################################################################################################################################################
"""

    flow_out(lf::Union{Leaf{FT}, Leaves2D{FT}}) where {FT}

Return the net flow that escape from the leaf, given
- `lf` `Leaf`, `Leaves2D`, `Root`, or `Stem` type organ

"""
function flow_out end

flow_out(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT} = flow_out(organ.HS.FLOW);

flow_out(mode::SteadyStateFlow{FT}) where {FT} = mode.flow;

flow_out(mode::NonSteadyStateFlow{FT,DIM_XYLEM}) where {FT,DIM_XYLEM} = mode.f_out;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to new version
#     2022-May-27: add method for steady flow mode
#     2022-May-27: add method for non-steady flow mode
#     2022-May-27: add method for root hydraulic system
#     2022-Jul-08: add method for root organ
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2023-Mar-28: return (NaN,0) if root is disconnected
#
#######################################################################################################################################################################################################
"""

    root_pk(root::Root{FT}) where {FT}

Return the root end pressure and total hydraulic conductance to find solution of flow rates in all roots, given
- `root` `Root` type struct

"""
function root_pk end

root_pk(root::Root{FT}, slayer::SoilLayer{FT}) where {FT} = root._isconnected ? root_pk(root.HS, slayer, root.t) : (FT(NaN), FT(0));

root_pk(hs::RootHydraulics{FT,DIM_XYLEM}, slayer::SoilLayer{FT}, T::FT) where {FT,DIM_XYLEM} = root_pk(hs, slayer, hs.FLOW, T);

root_pk(hs::RootHydraulics{FT,DIM_XYLEM}, slayer::SoilLayer{FT}, mode::SteadyStateFlow{FT}, T::FT) where {FT,DIM_XYLEM} = (
    (; AREA, K_RHIZ, K_X, L, ΔH) = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;
    _r_all::FT = 0;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _ in 1:10
        _k = relative_hydraulic_conductance(slayer.VC, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.flow / _k;
        _r_all += 1 / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T₂₅(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * _k_max * DIM_XYLEM;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= mode.flow / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;
        _r_all += 1 / _k;
    end;

    return _p_end, 1/_r_all
);

root_pk(hs::RootHydraulics{FT,DIM_XYLEM}, slayer::SoilLayer{FT}, mode::NonSteadyStateFlow{FT,DIM_XYLEM}, T::FT) where {FT,DIM_XYLEM} = (
    (; AREA, K_RHIZ, K_X, L, ΔH) = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;
    _r_all::FT = 0;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _ in 1:10
        _k = relative_hydraulic_conductance(slayer.VC, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.f_in / _k;
        _r_all += 1 / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T₂₅(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * _k_max * DIM_XYLEM;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= mode._f_element[_i] / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;
        _r_all += 1 / _k;
    end;

    return _p_end, 1/_r_all
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to new version
#     2022-May-27: rename the functions (flow_profile! and update_PVF!) to xylem_flow_profile!
#     2022-May-27: add method to set up root flow rate at steady state mode (this function is used to solve for steady state solution)
#     2022-May-27: add method to set up root flow rate at non-steady state mode (this function is used to solve for steady state solution)
#     2022-May-27: add method for leaf, root, and stem at steady state mode
#     2022-May-27: add method for leaf at non-steady state mode
#     2022-May-27: add method for root and stem at non-steady state mode
#     2022-May-27: add method for leaf, root, and stem hydraulic system at steady and non-steady state mode (for dispatching purpose)
#     2022-May-27: add method for leaf, root, and stem organ at steady and non-steady state mode (for dispatching purpose)
#     2022-May-27: add method to solve root flow rate partition at both steady and non-steady state modes
#     2022-May-27: add method for MonoElementSPAC (blank)
#     2022-May-31: remove hydraulic system from input variables, thus supporting leaf and stem
#     2022-May-31: use reformulate methods for setting up flow rate
#     2022-May-31: set up the flow rate profile using the network
#     2022-May-31: add method for MonoGrassSPAC
#     2022-May-31: add method for MonoPalmSPAC
#     2022-May-31: add method for MonoTreeSPAC
#     2022-Jun-29: rename SPAC to ML*SPAC to be more accurate
#     2022-Jun-30: add support to Leaves2D
#     2022-Jun-30: add method for Leaves1D
#     2022-Jul-12: add method to update leaf hydraulic flow rates per canopy layer based on stomatal conductance
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2022-Oct-20: fix a bug in flow profile counter (does not impact simulation)
#     2022-Oct-21: add a second solver to fix the case when root_pk does not work
#     2023-Mar-28: if root is disconnected, do not update its flow profile
#     2023-Mar-28: if no root is connected to soil, set all flow to 0
#
#######################################################################################################################################################################################################
"""
This function is designed to serve the following functionalities:
- Update flow profile in different organs
- Partition root flow rates at different layers
- Update flow profile for entire SPAC

"""
function xylem_flow_profile! end


"""

    xylem_flow_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Stem{FT}}, Δt::FT) where {FT}
    xylem_flow_profile!(organ::Root{FT}, Δt::FT) where {FT}
    xylem_flow_profile!(organ::Leaves1D{FT}, Δt::FT) where {FT}

Update organ flow rate profile after setting up the flow rate out, given
- `organ` `Leaf`, `Leaves1D`, `Leaves2D`, `Root`, or `Stem` type struct
- `Δt` Time step length

"""
xylem_flow_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Stem{FT}}, Δt::FT) where {FT} = xylem_flow_profile!(organ.HS, organ.t, Δt);

xylem_flow_profile!(organ::Root{FT}, Δt::FT) where {FT} = organ._isconnected ? xylem_flow_profile!(organ.HS, organ.t, Δt) : nothing;

xylem_flow_profile!(organ::Leaves1D{FT}, Δt::FT) where {FT} = (
    (; HS, HS2) = organ;

    xylem_flow_profile!(HS, organ.t[1], Δt);
    xylem_flow_profile!(HS2, organ.t[2], Δt);

    return nothing
);

xylem_flow_profile!(hs::Union{LeafHydraulics{FT,DIM_XYLEM}, RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}, T::FT, Δt::FT) where {FT,DIM_XYLEM} = xylem_flow_profile!(hs, hs.FLOW, T, Δt);

xylem_flow_profile!(hs::Union{LeafHydraulics{FT,DIM_XYLEM}, RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}, mode::SteadyStateFlow{FT}, T::FT, Δt::FT) where {FT,DIM_XYLEM} = nothing;

xylem_flow_profile!(hs::LeafHydraulics{FT,DIM_XYLEM}, mode::NonSteadyStateFlow{FT,1}, T::FT, Δt::FT) where {FT,DIM_XYLEM} = (
    (; PVC, V_MAXIMUM) = hs;

    _f_vis = relative_viscosity(T);

    # compute the flow rate from capacitance buffer
    mode._f_buffer[1] = (hs._p_storage - hs.p_leaf) * capacitance_buffer(PVC) / _f_vis * V_MAXIMUM;

    # make sure the buffer rate does not drain or overflow the capacictance
    if (mode._f_buffer[1] > 0) && (hs.v_storage <= mode._f_buffer[1] * Δt)
        mode._f_buffer[1] = (hs.v_storage - eps(FT)) / Δt;
    end;

    # update storage and the tissue pressure (p_storage)
    hs.v_storage -= mode._f_buffer[1] * Δt;
    hs._p_storage = xylem_pressure(PVC, hs.v_storage/V_MAXIMUM, T);

    # update flow into the tissue
    mode.f_in = mode.f_out - mode._f_buffer[1];

    return nothing
);

xylem_flow_profile!(hs::Union{RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}, mode::NonSteadyStateFlow{FT,DIM_XYLEM}, T::FT, Δt::FT) where {FT,DIM_XYLEM} = (
    (; PVC, V_MAXIMUM) = hs;

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


"""

    xylem_flow_profile!(roots::Vector{Root{FT}}, soil::Soil{FT}, cache_f::Vector{FT}, cache_k::Vector{FT}, cache_p::Vector{FT}, f_sum::FT, Δt::FT) where {FT}

Partition root flow rates at different layers for known total flow rate out, given
- `roots` Vector of `Root` in a multiple roots system
- `roots_index` Vector to match roots to soil layers
- `soil` Soil of companion roots
- `cache_f` Flow rate cache into each root
- `cache_k` Total conductance cache of each root
- `cache_p` Root xylem end pressure cache of each root
- `f_sum` Total flow rate out of the roots
- `Δt` Time step length

"""
xylem_flow_profile!(spac::MultiLayerSPAC{FT}, f_sum::FT, Δt::FT) where {FT} = (
    (; BRANCHES, LEAVES, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    # very first step here: if soil is too dry, disconnect root from soil
    _connected = 0;
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        _ψ_soil = soil_ψ_25(_slayer.VC, _slayer.θ) * relative_surface_tension(_slayer.t);
        _p_crit = critical_pressure(_root.HS.VC) * relative_surface_tension(_root.t);
        if _ψ_soil <= _p_crit
            disconnect!(_root);
        else
            _root._isconnected = true;
            _connected += 1;
        end;
    end;

    # if all roots are disconnected, set all flows to 0
    if _connected > 0
        spac._root_connection = true;
    else
        spac._root_connection = false;
        disconnect!(TRUNK);
        disconnect!.(BRANCHES);
        disconnect!.(LEAVES);

        return nothing
    end;

    # update root buffer rates to get an initial guess (flow rate not changing now as time step is set to 0)
    xylem_flow_profile!.(ROOTS, FT(0));

    # recalculate the flow profiles to make sure sum are the same as f_sum
    _use_second_solver = false;
    for _count in 1:20
        # sync the values to ks, ps, and qs
        for _i in eachindex(ROOTS_INDEX)
            _root = ROOTS[_i];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
            if _root._isconnected
                xylem_flow_profile!(_root.HS.FLOW, spac._fs[_i]);
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
        _f_diff = sum(spac._fs) - f_sum;
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

    # use second solver to solve for the flow rates (when SWC differs alot among layers)
    if _use_second_solver
        @inline diff_p_root(ind::Int, e::FT, p::FT) where {FT} = (
            _root = ROOTS[ind];
            _slayer = SOIL.LAYERS[ROOTS_INDEX[ind]];
            if _root._isconnected
                xylem_flow_profile!(ROOTS[ind].HS.FLOW, e);
            end;
            (_p,_) = root_pk(_root, _slayer);

            return _p - p
        );

        @inline diff_e_root(p::FT) where {FT} = (
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

            return _sum - f_sum
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
                xylem_flow_profile!(ROOTS[_i].HS.FLOW, spac._fs[_i]);
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
            xylem_flow_profile!(ROOTS[_i].HS.FLOW, spac._fs[_i]);
        end;
    end;
    xylem_flow_profile!.(ROOTS, Δt);

    return nothing
);

xylem_flow_profile!(mode::SteadyStateFlow, flow::FT) where {FT} = (mode.flow = flow; return nothing);

xylem_flow_profile!(mode::NonSteadyStateFlow, f_out::FT) where {FT} = (mode.f_out = f_out; mode._f_element .= f_out .- mode._f_sum; return nothing);


"""

    xylem_flow_profile!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT}
    xylem_flow_profile!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT}

Update flow profiles for the soil-plant-air continuum (set up leaf flow rate from stomatal conductance first), given
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type SPAC system
- `Δt` Time step length

"""
xylem_flow_profile!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT} = (
    (; LEAF, ROOT, STEM) = spac;

    # 0. update leaf flow or f_out from stomatal conductance
    xylem_flow_profile!(spac);

    # 1. update the leaf flow profile
    xylem_flow_profile!(LEAF, Δt);

    # 2. set up stem flow rate and profile
    xylem_flow_profile!(STEM.HS.FLOW, flow_in(LEAF) * LEAF.HS.AREA);
    xylem_flow_profile!(STEM, Δt);

    # 3. set up root flow rate and profile
    xylem_flow_profile!(ROOT.HS.FLOW, flow_in(STEM));
    xylem_flow_profile!(ROOT, Δt);

    return nothing
);

xylem_flow_profile!(spac::MultiLayerSPAC{FT}, Δt::FT) where {FT} = (
    (; BRANCHES, LEAVES, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    # 0. update leaf flow or f_out from stomatal conductance
    xylem_flow_profile!(spac);

    # 1. update the leaf flow profile
    xylem_flow_profile!.(LEAVES, Δt);

    # 2. set up branch flow rate and profile
    for _i in eachindex(LEAVES)
        xylem_flow_profile!(BRANCHES[_i].HS.FLOW, flow_in(LEAVES[_i]) * LEAVES[_i].HS.AREA);
    end;
    xylem_flow_profile!.(BRANCHES, Δt);

    # 3. set up trunk flow rate and profile
    xylem_flow_profile!(TRUNK.HS.FLOW, flow_in(BRANCHES));
    xylem_flow_profile!(TRUNK, Δt);

    # 4. set up root flow rate and profile
    xylem_flow_profile!(spac, flow_in(TRUNK), Δt);

    return nothing
);

xylem_flow_profile!(spac::MonoElementSPAC{FT}) where {FT} = (
    (; AIR, LEAF) = spac;

    # update the
    _g = 1 / (1 / LEAF.g_H₂O_s + 1 / (FT(1.35) * LEAF.g_CO₂_b));
    _d = saturation_vapor_pressure(LEAF.t) - AIR.p_H₂O;
    _f = _g * _d / AIR.P_AIR;
    xylem_flow_profile!(LEAF.HS.FLOW, _f);

    return nothing
);

xylem_flow_profile!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; AIR, CANOPY, DIM_LAYER, LEAVES, LEAVES_INDEX) = spac;

    for _i in eachindex(LEAVES)
        _p_sl = CANOPY.OPTICS.p_sunlit[DIM_LAYER + 1 - _i];

        _g_sh = 1 / (1 /LEAVES[_i].g_H₂O_s_shaded + 1 / (FT(1.35) * LEAVES[_i].g_CO₂_b));
        _g_sl = 0;
        for _j in eachindex(LEAVES[_i].g_H₂O_s_sunlit)
            _g_sl += 1 / (1 /LEAVES[_i].g_H₂O_s_sunlit[_j] + 1 / (FT(1.35) * LEAVES[_i].g_CO₂_b));
        end;
        _g_sl /= length(LEAVES[_i].g_H₂O_s_sunlit);

        _g = _g_sh * (1 - _p_sl) + _g_sl * _p_sl;
        _d = saturation_vapor_pressure(LEAVES[_i].t) - AIR[LEAVES_INDEX[_i]].p_H₂O;
        _f = _g * _d / AIR[LEAVES_INDEX[_i]].P_AIR;

        xylem_flow_profile!(LEAVES[_i].HS.FLOW, _f);
    end;

    return nothing
);

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
#     2023-Sep-07: update flow rate integrators
#     2023-Sep-07: add root disconnection warning message
#
#######################################################################################################################################################################################################
"""

    root_pk(root::Root{FT}) where {FT<:AbstractFloat}

Return the root end pressure and total hydraulic conductance to find solution of flow rates in all roots, given
- `root` `Root` type struct

"""
function root_pk end

root_pk(root::Root{FT}, slayer::SoilLayer{FT}) where {FT<:AbstractFloat} = root._isconnected ? root_pk(root.HS, slayer, root.t) : (FT(NaN), FT(0));

root_pk(hs::RootHydraulics{FT}, slayer::SoilLayer{FT}, T::FT) where {FT<:AbstractFloat} = root_pk(hs, slayer, hs.FLOW, T);

root_pk(hs::RootHydraulics{FT}, slayer::SoilLayer{FT}, mode::SteadyStateFlow{FT}, T::FT) where {FT<:AbstractFloat} = (
    (; AREA, DIM_XYLEM, K_RHIZ, K_X, L, ΔH) = hs;

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

root_pk(hs::RootHydraulics{FT}, slayer::SoilLayer{FT}, mode::NonSteadyStateFlow{FT}, T::FT) where {FT<:AbstractFloat} = (
    (; AREA, DIM_XYLEM, K_RHIZ, K_X, L, ΔH) = hs;

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

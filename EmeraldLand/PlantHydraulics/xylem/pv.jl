# This file contains the pressure volume curve functions

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-24: add method for LinearPVCurve
#     2022-May-24: add method for SegmentedPVCurve
#     2023-Sep-22: add method for ExponentialPVCurve
#
#######################################################################################################################################################################################################
"""

    capacitance_buffer(pv::ExponentialPVCurve{FT}, rvol::FT, t::FT) where {FT}
    capacitance_pressure(pv::LinearPVCurve{FT}, rvol::FT, t::FT) where {FT}
    capacitance_pressure(pv::SegmentedPVCurve{FT}, rvol::FT, t::FT) where {FT}

Return the xylem water pressure in MPa, given
- `pv` `ExponentialPVCurve`, `LinearPVCurve` or `SegmentedPVCurve` type pressure volume curve
- `rvol` Relative volume (relative to maximum)
- `t` Temperature

"""
function capacitance_pressure end;

capacitance_pressure(pv::ExponentialPVCurve{FT}, rvol::FT, t::FT) where {FT} = log(rvol) / pv.slope;

capacitance_pressure(pv::LinearPVCurve{FT}, rvol::FT, t::FT) where {FT} = (rvol - 1) / pv.slope;

capacitance_pressure(pv::SegmentedPVCurve{FT}, rvol::FT, t::FT) where {FT} = (
    if rvol > pv.θ_tlp
        return -pv.c_all * GAS_R(FT) * t / (rvol - pv.θ_apo) * FT(1e-6) + pv.ϵ_bulk * (rvol - pv.θ_tlp)
    elseif rvol > pv.θ_apo
        return max(-100, -pv.c_all * GAS_R(FT) * t / (rvol - pv.θ_apo) * FT(1e-6))
    else
        return FT(-100)
    end;
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-23: add function to compute the PV volume from pressure
#
#######################################################################################################################################################################################################
"""

    capacitance_volume(pv::ExponentialPVCurve{FT}, p::FT, t::FT) where {FT}
    capacitance_volume(pv::LinearPVCurve{FT}, p::FT, t::FT) where {FT}
    capacitance_volume(pv::SegmentedPVCurve{FT}, p::FT, t::FT) where {FT}

Return the relative capaciatance volume, given
- `pv` `ExponentialPVCurve`, `LinearPVCurve` or `SegmentedPVCurve` type pressure volume curve
- `p` Capacitance water pressure in MPa
- `t` Temperature

"""
function capacitance_volume end;

# TODO: if pressure is really low, the volume gets too close to zero, thus there is a problem. Consider adding a baseline here, say 0.2 + 0.8 * exp(p * pv.slope)
capacitance_volume(pv::ExponentialPVCurve{FT}, p::FT, t::FT) where {FT} = exp(p * pv.slope);

capacitance_volume(pv::LinearPVCurve{FT}, p::FT, t::FT) where {FT} = max(0, 1 + p * pv.slope);

capacitance_volume(pv::SegmentedPVCurve{FT}, p::FT, t::FT) where {FT} = (
    p_tlp = -pv.c_all * GAS_R(FT) * t / (pv.θ_tlp - pv.θ_apo) * FT(1e-6);
    if p <= p_tlp
        return pv.θ_apo - pv.c_all * GAS_R(FT) * t * FT(1e-6) / p
    else
        a = pv.ϵ_bulk;
        b = -1 * (pv.ϵ_bulk * (pv.θ_tlp + pv.θ_apo) + p);
        c = pv.ϵ_bulk * pv.θ_tlp * pv.θ_apo + p * pv.θ_apo - pv.c_all * GAS_R(FT) * t * FT(1e-6);

        return upper_quadratic(a, b, c)
    end;
);

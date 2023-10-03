
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function for StomataModels.jl (for nocturnal transpiration)
#     2022-Jul-29: add support to Q10Peak
#
#######################################################################################################################################################################################################
"""

    ∂R∂T(leaf::Leaf{FT}) where {FT}
    ∂R∂T(leaves::Leaves2D{FT}) where {FT}

Return the marginal increase in respiration rate per temperature, given
- `leaf` `Leaf` type leaf
- `leaves` `Leaves2D` type leaf

"""
function ∂R∂T end

∂R∂T(leaves::Leaves2D{FT}) where {FT} = ∂R∂T(leaves.PSM, leaves.t);

∂R∂T(psm::Union{C3Cyto{FT}, C3VJP{FT}, C4VJP{FT}}, t::FT) where {FT} = ∂R∂T(psm.TD_R, psm.r_d25, t);

∂R∂T(td::Arrhenius{FT}, r_ref::FT, t::FT) where {FT} = r_ref * exp(td.ΔHA / GAS_R(FT) * (1/td.T_REF - 1/t)) * td.ΔHA / (GAS_R(FT) * t ^ 2);

∂R∂T(td::ArrheniusPeak{FT}, r_ref::FT, t::FT) where {FT} = (
    (; T_REF, ΔHA, ΔHD, ΔSV) = td;

    # _f_a: activation correction, _f_b: de-activation correction
    _expt = exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t));
    _rt²  = GAS_R(FT) * t ^ 2;
    _f_a  = exp(ΔHA / GAS_R(FT) * (1 / T_REF - 1 / t));
    _f_b  = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * T_REF))) / (1 + _expt);
    _f_a′ = _f_a * ΔHA / _rt²;
    _f_b′ = -1 * _f_b / (1 + _expt) * _expt * ΔHD / _rt²;

    return r_ref * (_f_a′ * _f_b + _f_a * _f_b′)
);

∂R∂T(td::Q10{FT}, r_ref::FT, t::FT) where {FT} = r_ref * log(td.Q_10) * td.Q_10 ^ ( (t - td.T_REF) / 10) / 10;

∂R∂T(td::Q10Peak{FT}, r_ref::FT, t::FT) where {FT} = (
    (; T_REF, ΔHD, ΔSV) = td;

    # _f_a: activation correction, _f_b: de-activation correction
    _expt = exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t));
    _rt²  = GAS_R(FT) * t ^ 2;
    _f_a  = td.Q_10 ^ ( (t - T_REF) / 10 );
    _f_b  = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * T_REF))) / (1 + _expt);
    _f_a′ = log(td.Q_10) * td.Q_10 ^ ( (t - td.T_REF) / 10) / 10;
    _f_b′ = -1 * _f_b / (1 + _expt) * _expt * ΔHD / _rt²;

    return r_ref * (_f_a′ * _f_b + _f_a * _f_b′)
);

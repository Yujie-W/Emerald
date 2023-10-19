# This file contains functions to compute the derivatives to use with StomataModels.jl

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

Return the marginal increase in respiration rate per temperature, given
- `leaf` `Leaf` type leaf

"""
function ∂R∂T end;

∂R∂T(leaf::Leaf{FT}) where {FT} = ∂R∂T(leaf.photosystem, leaf.energy.auxil.t);

∂R∂T(ps::Union{C3Cyto{FT}, C3VJP{FT}, C4VJP{FT}}, t::FT) where {FT} = ∂R∂T(ps.state.TD_R, ps.state.r_d25, t);

∂R∂T(td::Arrhenius{FT}, r_ref::FT, t::FT) where {FT} = r_ref * exp(td.ΔHA / GAS_R(FT) * (1/td.T_REF - 1/t)) * td.ΔHA / (GAS_R(FT) * t ^ 2);

∂R∂T(td::ArrheniusPeak{FT}, r_ref::FT, t::FT) where {FT} = (
    (; T_REF, ΔHA, ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    expt = exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t));
    rt²  = GAS_R(FT) * t ^ 2;
    f_a  = exp(ΔHA / GAS_R(FT) * (1 / T_REF - 1 / t));
    f_b  = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * T_REF))) / (1 + expt);
    f_a′ = f_a * ΔHA / rt²;
    f_b′ = -1 * f_b / (1 + expt) * expt * ΔHD / rt²;

    return r_ref * (f_a′ * f_b + f_a * f_b′)
);

∂R∂T(td::Q10{FT}, r_ref::FT, t::FT) where {FT} = r_ref * log(td.Q_10) * td.Q_10 ^ ( (t - td.T_REF) / 10) / 10;

∂R∂T(td::Q10Peak{FT}, r_ref::FT, t::FT) where {FT} = (
    (; T_REF, ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    expt = exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t));
    rt²  = GAS_R(FT) * t ^ 2;
    f_a  = td.Q_10 ^ ( (t - T_REF) / 10 );
    f_b  = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * T_REF))) / (1 + expt);
    f_a′ = log(td.Q_10) * td.Q_10 ^ ( (t - td.T_REF) / 10) / 10;
    f_b′ = -1 * f_b / (1 + expt) * expt * ΔHD / rt²;

    return r_ref * (f_a′ * f_b + f_a * f_b′)
);

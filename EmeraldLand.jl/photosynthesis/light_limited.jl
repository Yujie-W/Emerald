#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function light_limited_rate!
#
#######################################################################################################################################################################################################
"""
This function supports two types of calculations:
- Calculate the rate from internal CO₂
- Calculate the rate from CO₂ conductance by solving a quadratic function

"""
function light_limited_rate! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: refactor the function light_limited_rate!
#     2022-Jan-14: add p_i to input list to make the code more modular
#     2022-Jan-24: add C3CytochromeModel support in a Union
#     2022-Feb-07: remove duplicated j (using j_pot is enough)
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}
    light_limited_rate!(psm::C3VJPModel{FT}) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` structure for C3 photosynthesis model

"""
light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat} = (psm._a_j = psm._j_pot * psm._e_to_c; return nothing);

light_limited_rate!(psm::C3VJPModel{FT}) where {FT<:AbstractFloat} = (psm._a_j = psm._j * psm._e_to_c; return nothing);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add g_lc to input list to make the code more modular
#     2022-Jan-14: add this new method to simplify the multiple dispatch of leaf_photosynthesis!
#     2022-Feb-07: remove duplicated j (using j_pot is enough for C4VJPModel)
#     2022-Feb-28: add C3CytochromeModel support
#     2022-Mar-01: j_psi, η_c, and η_l from psm (temperature corrected) rather than constant Η_C and Η_L
#     2022-Jun-27: remove ppar from input variable list
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    light_limited_rate!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    light_limited_rate!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` structure for C3 photosynthesis model
- `rc` `CytochromeReactionCenter` or `VJPReactionCenter` type photosynthesis system reaction center
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    (; EFF_1, EFF_2) = psm;

    _eff_a = 1 - psm._η_l / psm._η_c;
    _eff_b = 1 / psm._η_c;
    _eff_1 = _eff_a * EFF_1 + 3 * _eff_b;
    _eff_2 = _eff_a * EFF_2 + 7 * _eff_b;

    _a = psm._j_psi;
    _b = psm._j_psi * psm._γ_star;
    _c = _eff_1;
    _d = _eff_2 * psm._γ_star;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = β * psm._r_d;

    _qa = _c * _f;
    _qb = _c*_f*_r - _c*_p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_c*_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    psm._a_j = _an + _r;

    return nothing
);

light_limited_rate!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _a = psm._j;
    _b = psm._j * psm._γ_star;
    _c = psm.EFF_1;
    _d = psm.EFF_2 * psm._γ_star;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = β * psm._r_d;

    _qa = _c * _f;
    _qb = _c*_f*_r - _c*_p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_c*_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    psm._a_j = _an + _r;

    return nothing
);

light_limited_rate!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (psm._a_j = psm._j_pot * psm._e_to_c; return nothing);

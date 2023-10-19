#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function product_limited_rate!
#     2022-Jan-14: add input variable p_i to make the code more modular
#     2022-Jan-14: add input variable g_lc to make the code more modular
#     2022-Jan-18: add support to C3CytochromeModel
#     2022-Feb-28: add C3CytochromeModel support
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
# Bug fixes
#     2023-Sep-21: if g_lc is 0, set a_p to r
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, p_i::FT; β::FT = FT(1)) where {FT}
    product_limited_rate!(psm::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT}

Update the product limited photosynthetic rate (p_i for PCO₂Mode, g_lc for GCO₂Mode), given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `air` `AirLayer` struct
- `g_lc` Leaf conductance in `mol m⁻² s⁻¹`

"""
function product_limited_rate! end;

product_limited_rate!(psm::Union{C3Cyto{FT}, C3VJP{FT}}, p_i::FT; β::FT = FT(1)) where {FT} = (psm.auxil.a_p = β * psm.auxil.v_cmax / 2; return nothing);

product_limited_rate!(psm::C4VJP{FT}, p_i::FT; β::FT = FT(1)) where {FT} = (psm.auxil.a_p = β * psm.auxil.v_pmax * p_i / (p_i + psm.auxil.k_pep); return nothing);

product_limited_rate!(psm::Union{C3Cyto{FT}, C3VJP{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (psm.auxil.a_p = β * psm.auxil.v_cmax / 2; return nothing);

product_limited_rate!(psm::C4VJP{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    a = β * psm.auxil.v_pmax;
    d = psm.auxil.k_pep;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.auxil.ps[2];
    r = β * psm.auxil.r_d;

    qa = f;
    qb = f * r - p - d - a * f;
    qc = a * p - r * (p + d);
    an = lower_quadratic(qa, qb, qc);

    if g_lc == 0
        psm.auxil.a_p = r;
    else
        psm.auxil.a_p = an + r;
    end;

    return nothing
);

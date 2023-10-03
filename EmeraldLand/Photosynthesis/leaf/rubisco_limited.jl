# This file contains functions to compute RubisCO limited photosynthesis rate

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function rubisco_limited_rate!
#     2022-Jan-14: add input variable p_i to make the code more modular
#     2022-Jan-14: add input variable g_lc to make the code more modular
#     2022-Feb-07: add support to C3CytochromeModel
#     2022-Feb-28: add C3CytochromeModel support
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
# Bug fixes
#     2023-Sep-21: if g_lc is 0, set a_c to r
#
#######################################################################################################################################################################################################
"""

    rubisco_limited_rate!(psm::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, p_i::FT; β::FT = FT(1)) where {FT}
    rubisco_limited_rate!(psm::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT}

Update the RubisCO limited photosynthetic rate (p_i for PCO₂Mode and g_lc for GCO₂Mode), given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`

"""
function rubisco_limited_rate! end

rubisco_limited_rate!(psm::Union{C3Cyto{FT},C3VJP{FT}}, p_i::FT; β::FT = FT(1)) where {FT} = (
    psm.auxil.a_c = β * psm.auxil.v_cmax * (p_i - psm.auxil.γ_star) / (p_i + psm.auxil.k_m);

    return nothing
);

rubisco_limited_rate!(psm::C4VJP{FT}, p_i::FT; β::FT = FT(1)) where {FT} = (psm.auxil.a_c = β * psm.auxil.v_cmax; return nothing);

rubisco_limited_rate!(psm::Union{C3Cyto{FT}, C3VJP{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    a = β * psm.auxil.v_cmax;
    b = β * psm.auxil.v_cmax * psm.auxil.γ_star;
    d = psm.auxil.k_m;
    f = air.P_AIR / g_lc * FT(1e-6);
    p = air.p_CO₂;
    r = β * psm.auxil.r_d;

    qa = f;
    qb = f*r - p - d - a*f;
    qc = a*p - b - r*(p + d);
    an = lower_quadratic(qa, qb, qc);

    if g_lc == 0
        psm.auxil.a_c = r;
    else
        psm.auxil.a_c = an + r;
    end;

    return nothing
);

rubisco_limited_rate!(psm::C4VJP{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (psm.auxil.a_c = β * psm.auxil.v_cmax; return nothing);

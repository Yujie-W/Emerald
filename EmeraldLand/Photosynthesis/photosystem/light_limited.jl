#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function light_limited_rate!
#     2022-Jan-14: add p_i to input list to make the code more modular
#     2022-Jan-14: add g_lc to input list to make the code more modular
#     2022-Jan-24: add C3CytochromeModel support in a Union
#     2022-Feb-28: add C3CytochromeModel support
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Jun-15: set a_j to 0 when j is 0 (not a quadratic function any more)
# Bug fixes
#     2023-Sep-21: if g_lc is 0, set a_j to r
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::LeafPhotosystem{FT}) where {FT}
    light_limited_rate!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT}

Update the electron transport limited photosynthetic rate (none for PCO₂Mode and g_lc for GCO₂Mode), given
- `psm` `LeafPhotosystem` struct
- `air` `AirLayer` struct for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function light_limited_rate! end;

# For CanopyLayer
# PCO₂Mode
light_limited_rate!(psm::CanopyLayerPhotosystem{FT}) where {FT} = light_limited_rate!(psm.auxil);

light_limited_rate!(psa::CanopyLayerPhotosystemAuxil{FT}) where {FT} = (psa.a_j .= psa.j .* psa.e2c; return nothing);

# GCO₂Mode
light_limited_rate!(psm::CanopyLayerPhotosystem{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} = light_limited_rate!(psm.trait, psm.state, psm.auxil, air, g_lc; β = β);

light_limited_rate!(
            pst::Union{C3CytoMinEtaTrait{FT}, C3CytoTrait{FT}, C3JBTrait{FT}},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if psa.j_psi[1] == 0 && psa.j_psi[end] == 0
        psa.a_j .= 0;

        return nothing
    end;

    eff_a = 1 - psa.η_l / psa.η_c;
    eff_b = 1 / psa.η_c;
    eff_1 = eff_a * pss.EFF_1 + 3 * eff_b;
    eff_2 = eff_a * pss.EFF_2 + 7 * eff_b;

    a = psa.j_psi;
    b = psa.j_psi * psa.γ_star;
    c = eff_1;
    d = eff_2 * psa.γ_star;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = c * f;
    qb = c * f * r - c * p - d - a * f;
    qc = a * p - b - r * (c * p + d);
    an = lower_quadratic.(qa, qb, qc);

    if g_lc[1] == 0 && g_lc[end] == 0
        psa.a_j .= r;
    else
        psa.a_j .= an .+ r;
    end;

    return nothing
);

light_limited_rate!(
            pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if psa.j[1] == 0 && psa.j[end] == 0
        psa.a_j .= 0;

        return nothing
    end;

    a = psa.j;
    b = psa.j * psa.γ_star;
    c = pss.EFF_1;
    d = pss.EFF_2 * psa.γ_star;
    f = air.state.p_air ./ g_lc .* FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = c .* f;
    qb = c .* f .* r .- c .* p .- d .- a .* f;
    qc = a .* p .-  b .- r .* (c .* p .+ d);
    an = lower_quadratic.(qa, qb, qc);

    if g_lc[1] == 0 && g_lc[end] == 0
        psa.a_j .= r;
    else
        psa.a_j .= an .+ r;
    end;

    return nothing
);

light_limited_rate!(pst::Union{C4CLMTrait{FT}, C4VJPTrait{FT}}, pss::C4State{FT}, psa::CanopyLayerPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} =
    (psa.a_j .= psa.j .* psa.e2c; return nothing);


# For Leaf
# PCO₂Mode
light_limited_rate!(psm::LeafPhotosystem{FT}) where {FT} = light_limited_rate!(psm.auxil);

light_limited_rate!(psa::LeafPhotosystemAuxil{FT}) where {FT} = (psa.a_j = psa.j * psa.e2c; return nothing);

# GCO₂Mode
light_limited_rate!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = light_limited_rate!(psm.trait, psm.state, psm.auxil, air, g_lc; β = β);

light_limited_rate!(pst::Union{C3CytoMinEtaTrait{FT}, C3CytoTrait{FT}, C3JBTrait{FT}}, pss::C3State{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    if psa.j_psi == 0
        psa.a_j = 0;

        return nothing
    end;

    eff_a = 1 - psa.η_l / psa.η_c;
    eff_b = 1 / psa.η_c;
    eff_1 = eff_a * pss.EFF_1 + 3 * eff_b;
    eff_2 = eff_a * pss.EFF_2 + 7 * eff_b;

    a = psa.j_psi;
    b = psa.j_psi * psa.γ_star;
    c = eff_1;
    d = eff_2 * psa.γ_star;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = c * f;
    qb = c * f * r - c * p - d - a * f;
    qc = a * p - b - r * (c * p + d);
    an = lower_quadratic(qa, qb, qc);

    if g_lc == 0
        psa.a_j = r;
    else
        psa.a_j = an + r;
    end;

    return nothing
);

light_limited_rate!(pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}}, pss::C3State{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    if psa.j == 0
        psa.a_j = 0;

        return nothing
    end;

    a = psa.j;
    b = psa.j * psa.γ_star;
    c = pss.EFF_1;
    d = pss.EFF_2 * psa.γ_star;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = c * f;
    qb = c * f * r - c * p - d - a * f;
    qc = a * p -  b - r * (c * p + d);
    an = lower_quadratic(qa, qb, qc);

    if g_lc == 0
        psa.a_j = r;
    else
        psa.a_j = an + r;
    end;

    return nothing
);

light_limited_rate!(pst::Union{C4CLMTrait{FT}, C4VJPTrait{FT}}, pss::C4State{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} =
    (psa.a_j = psa.j * psa.e2c; return nothing);

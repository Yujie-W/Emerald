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
#     2024-Apr-15: add support to C4CLMTrait model using CLM settings
#     2024-Jul-22: add support to C3FvCBTrait and C3JBTrait model (infinity a_p)
#     2024-Jul-23: add support to C3CytoInfApTrait model (infinity a_p)
# Bug fixes
#     2023-Sep-21: if g_lc is 0, set a_p to r
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::LeafPhotosystem{FT}, p_i::FT; β::FT = FT(1)) where {FT}
    product_limited_rate!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT}

Update the product limited photosynthetic rate (p_i for PCO₂Mode, g_lc for GCO₂Mode), given
- `psm` `LeafPhotosystem` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `air` `AirLayer` struct
- `g_lc` Leaf conductance in `mol m⁻² s⁻¹`

"""
function product_limited_rate! end;

# For CanopyLayer
# PCO₂Mode
product_limited_rate!(psm::CanopyLayerPhotosystem{FT}, p_i::Vector{FT}; β::FT = FT(1)) where {FT} = product_limited_rate!(psm.trait, psm.auxil, p_i; β = β);

product_limited_rate!(pst::Union{C3CLMTrait{FT}, C3CytoTrait{FT}, C3VJPTrait{FT}}, psa::CanopyLayerPhotosystemAuxil{FT}, p_i::Vector{FT}; β::FT = FT(1)) where {FT} = (
    psa.a_p .= β .* psa.v_cmax ./ 2;

    return nothing
);

product_limited_rate!(pst::Union{C3CytoInfApTrait{FT}, C3FvCBTrait{FT}, C3JBTrait{FT}}, psa::CanopyLayerPhotosystemAuxil{FT}, p_i::Vector{FT}; β::FT = FT(1)) where {FT} = (
    psa.a_p .= FT(Inf);

    return nothing
);

product_limited_rate!(pst::C4CLMTrait{FT}, psa::CanopyLayerPhotosystemAuxil{FT}, p_i::Vector{FT}; β::FT = FT(1)) where {FT} = (
    psa.a_p .= β .* psa.k_pep_clm .* pst.v_cmax25 .* p_i;

    return nothing
);

product_limited_rate!(pst::C4VJPTrait{FT}, psa::CanopyLayerPhotosystemAuxil{FT}, p_i::Vector{FT}; β::FT = FT(1)) where {FT} = (
    psa.a_p .= β .* psa.v_pmax .* p_i ./ (p_i .+ psa.k_pep);

    return nothing
);

# GCO₂Mode
product_limited_rate!(psm::CanopyLayerPhotosystem{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} = product_limited_rate!(psm.trait, psm.auxil, air, g_lc; β = β);

product_limited_rate!(pst::Union{C3CLMTrait{FT}, C3CytoTrait{FT}, C3VJPTrait{FT}}, psa::CanopyLayerPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} = (
    psa.a_p .= β .* psa.v_cmax ./ 2;

    return nothing
);

product_limited_rate!(pst::Union{C3CytoInfApTrait{FT}, C3FvCBTrait{FT}, C3JBTrait{FT}}, psa::CanopyLayerPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} = (
    psa.a_p .= FT(Inf);

    return nothing
);

product_limited_rate!(pst::C4CLMTrait{FT}, psa::CanopyLayerPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} = (
    a = air.state.p_air;
    g = FT(1e6) * g_lc;
    k = β * psa.k_pep_clm * pst.v_cmax25;
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    p_i = (g * p + a * r) / (a * k + g);
    psa.a_p .= k .* p_i;

    return nothing
);

product_limited_rate!(pst::C4VJPTrait{FT}, psa::CanopyLayerPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::Vector{FT}; β::FT = FT(1)) where {FT} = (
    a = β * psa.v_pmax;
    d = psa.k_pep;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = f;
    qb = f * r - p - d - a * f;
    qc = a * p - r * (p + d);
    an = lower_quadratic.(qa, qb, qc);

    if g_lc[0] == 0 && g_lc[end] == 0
        psa.a_p .= r;
    else
        psa.a_p .= an .+ r;
    end;

    return nothing
);

# For Leaf
# PCO₂Mode
product_limited_rate!(psm::LeafPhotosystem{FT}, p_i::FT; β::FT = FT(1)) where {FT} = product_limited_rate!(psm.trait, psm.auxil, p_i; β = β);

product_limited_rate!(pst::Union{C3CLMTrait{FT}, C3CytoTrait{FT}, C3VJPTrait{FT}}, psa::LeafPhotosystemAuxil{FT}, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.a_p = β * psa.v_cmax / 2;

    return nothing
);

product_limited_rate!(pst::Union{C3CytoInfApTrait{FT}, C3FvCBTrait{FT}, C3JBTrait{FT}}, psa::LeafPhotosystemAuxil{FT}, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.a_p = FT(Inf);

    return nothing
);

product_limited_rate!(pst::C4CLMTrait{FT}, psa::LeafPhotosystemAuxil{FT}, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.a_p = β * psa.k_pep_clm * pst.v_cmax25 * p_i;

    return nothing
);

product_limited_rate!(pst::C4VJPTrait{FT}, psa::LeafPhotosystemAuxil{FT}, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.a_p = β * psa.v_pmax * p_i / (p_i + psa.k_pep);

    return nothing
);

# GCO₂Mode
product_limited_rate!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = product_limited_rate!(psm.trait, psm.auxil, air, g_lc; β = β);

product_limited_rate!(pst::Union{C3CLMTrait{FT}, C3CytoTrait{FT}, C3VJPTrait{FT}}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    psa.a_p = β * psa.v_cmax / 2;

    return nothing
);

product_limited_rate!(pst::Union{C3CytoInfApTrait{FT}, C3FvCBTrait{FT}, C3JBTrait{FT}}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    psa.a_p = FT(Inf);

    return nothing
);

product_limited_rate!(pst::C4CLMTrait{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    a = air.state.p_air;
    g = FT(1e6) * g_lc;
    k = β * psa.k_pep_clm * pst.v_cmax25;
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    p_i = (g * p + a * r) / (a * k + g);
    psa.a_p = k * p_i;

    return nothing
);

product_limited_rate!(pst::C4VJPTrait{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT} = (
    a = β * psa.v_pmax;
    d = psa.k_pep;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = f;
    qb = f * r - p - d - a * f;
    qc = a * p - r * (p + d);
    an = lower_quadratic(qa, qb, qc);

    if g_lc == 0
        psa.a_p = r;
    else
        psa.a_p = an + r;
    end;

    return nothing
);

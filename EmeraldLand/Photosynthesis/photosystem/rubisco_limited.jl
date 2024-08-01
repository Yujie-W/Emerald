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
#     2024-Aug-01: generalize the function for GeneralC3Trait and GeneralC4Trait
# Bug fixes
#     2023-Sep-21: if g_lc is 0, set a_c to r
#
#######################################################################################################################################################################################################
"""

    rubisco_limited_rate!(ps::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, p_i::FT; β::FT = FT(1)) where {FT}
    rubisco_limited_rate!(ps::Union{C3Cyto{FT},C3VJP{FT},C4VJP{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT}

Update the RubisCO limited photosynthetic rate, given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`

"""
function rubisco_limited_rate! end;

# For Leaf
# Pressure mode
rubisco_limited_rate!(
            ps::LeafPhotosystem{FT},
            p_i::FT;
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(ps.trait, ps.auxil, p_i; β = β);

rubisco_limited_rate!(
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            psa::LeafPhotosystemAuxil{FT},
            p_i::FT;
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(psa, pst.ACM, p_i; β = β);

rubisco_limited_rate!(
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC3VcmaxPi,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (psa.a_c = β * psa.v_cmax * (p_i - psa.γ_star) / (p_i + psa.k_m); return nothing);

rubisco_limited_rate!(
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC4Vcmax,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (psa.a_c = β * psa.v_cmax; return nothing);

# Conductance mode
rubisco_limited_rate!(
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(ps.state, ps.auxil, air, g_lc; β = β);

rubisco_limited_rate!(
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(psa, pst.ACM, air, g_lc; β = β);

rubisco_limited_rate!(
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC3VcmaxPi,
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = (
    a = β * psa.v_cmax;
    b = β * psa.v_cmax * psa.γ_star;
    d = psa.k_m;
    f = air.state.p_air / g_lc * FT(1e-6);
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    qa = f;
    qb = f*r - p - d - a*f;
    qc = a*p - b - r*(p + d);
    an = lower_quadratic(qa, qb, qc);

    if g_lc == 0
        psa.a_c = r;
    else
        psa.a_c = an + r;
    end;

    return nothing
);

rubisco_limited_rate!(
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC4Vcmax,
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = (psa.a_c = β * psa.v_cmax; return nothing);

# For CanopyLayer (Conductance mode only)
rubisco_limited_rate!(
            cache::SPACCache{FT},
            ps::CanopyLayerPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(cache, ps.trait, ps.auxil, air, g_lc; β = β);

rubisco_limited_rate!(
            cache::SPACCache{FT},
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(cache, psa, pst.ACM, air, g_lc; β = β);

rubisco_limited_rate!(
            cache::SPACCache{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            acm::AcMethodC3VcmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack the cache variables
    f  = cache.cache_incl_azi_2_1;
    qb = cache.cache_incl_azi_2_2;
    qc = cache.cache_incl_azi_2_3;
    an = cache.cache_incl_azi_2_4;

    a = β * psa.v_cmax;
    b = β * psa.v_cmax * psa.γ_star;
    d = psa.k_m;
    p = air.s_aux.ps[2];
    r = β * psa.r_d;
    @. f = air.state.p_air / g_lc * FT(1e-6);

    qa = f;
    @. qb = f * r - p - d - a * f;
    @. qc = a * p - b - r * (p + d);
    @. an = lower_quadratic(qa, qb, qc);

    if g_lc[1] == 0 && g_lc[end] == 0
        @. psa.a_c = r;
    else
        @. psa.a_c = an + r;
    end;

    return nothing
);

rubisco_limited_rate!(
            cache::SPACCache{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            acm::AcMethodC4Vcmax,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_c = β * psa.v_cmax; return nothing);

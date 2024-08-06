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
#     2024-Aug-01: generalize the function for GeneralC3Trait and GeneralC4Trait
# Bug fixes
#     2023-Sep-21: if g_lc is 0, set a_p to r
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(
                cache::SPACCache{FT},
                ps::CanopyLayerPhotosystem{FT},
                air::AirLayer{FT},
                g_lc::Vector{FT};
                β::FT = FT(1)) where {FT}
    product_limited_rate!(
                ps::LeafPhotosystem{FT},
                p_i::FT;
                β::FT = FT(1)) where {FT}
    product_limited_rate!(
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                g_lc::FT;
                β::FT = FT(1)) where {FT}

Update the product limited photosynthetic rate, given
- `cache` `SPACCache` struct
- `ps` `CanopyLayerPhotosystem` or `LeafPhotosystem` struct
- `air` `AirLayer` struct for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `p_i` Internal CO₂ partial pressure in `Pa`

"""
function product_limited_rate! end;

# For Leaf
# Pressure mode
product_limited_rate!(
            ps::LeafPhotosystem{FT},
            p_i::FT;
            β::FT = FT(1)) where {FT} = product_limited_rate!(ps.trait, ps.auxil, p_i; β = β);

product_limited_rate!(
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            psa::LeafPhotosystemAuxil{FT},
            p_i::FT;
            β::FT = FT(1)) where {FT} = product_limited_rate!(pst, psa, pst.APM, p_i; β = β);

product_limited_rate!(
            pst::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC3Inf,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (psa.a_p = FT(Inf); return nothing);

product_limited_rate!(
            pst::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC3Vcmax,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (psa.a_p = β * psa.v_cmax / 2; return nothing);

product_limited_rate!(
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC4VcmaxPi,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (psa.a_p = β * psa.k_pep_clm * pst.v_cmax25 * p_i; return nothing);

product_limited_rate!(
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC4VpmaxPi,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (psa.a_p = β * psa.v_pmax * p_i / (p_i + psa.k_pep); return nothing);

# Conductance mode
product_limited_rate!(
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = product_limited_rate!(ps.trait, ps.auxil, air, g_lc; β = β);

product_limited_rate!(
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = product_limited_rate!(pst, psa, pst.APM, air, g_lc; β = β);

product_limited_rate!(
            pst::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC3Inf,
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = (psa.a_p = FT(Inf); return nothing);

product_limited_rate!(
            pst::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ApMethodC3Vcmax,
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = (psa.a_p = β * psa.v_cmax / 2; return nothing);

product_limited_rate!(
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC4VcmaxPi,
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = (
    a = air.state.p_air;
    g = FT(1e6) * g_lc;
    k = β * psa.k_pep_clm * pst.v_cmax25;
    p = air.s_aux.ps[2];
    r = β * psa.r_d;

    p_i = (g * p + a * r) / (a * k + g);
    psa.a_p = k * p_i;

    return nothing
);

product_limited_rate!(
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            apm::ApMethodC4VpmaxPi,
            air::AirLayer{FT},
            g_lc::FT;
            β::FT = FT(1)) where {FT} = (
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

# For CanopyLayer (Conductance mode only)
product_limited_rate!(
            cache::SPACCache{FT},
            ps::CanopyLayerPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate!(cache, ps.trait, ps.auxil, air, g_lc; β = β);

product_limited_rate!(
            cache::SPACCache{FT},
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            air::AirLayer{FT}, g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate!(cache, pst, psa, pst.APM, air, g_lc; β = β);

product_limited_rate!(
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            apm::ApMethodC3Inf,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_p = FT(Inf); return nothing);

product_limited_rate!(
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            apm::ApMethodC3Vcmax,
            air::AirLayer{FT}, g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.v_cmax / 2; return nothing);

product_limited_rate!(
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            apm::ApMethodC4VcmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack variables from the cache
    g = cache.cache_incl_azi_2_1;
    p_i = cache.cache_incl_azi_2_2;

    a = air.state.p_air;
    k = β * psa.k_pep_clm * pst.v_cmax25;
    p = air.s_aux.ps[2];
    r = β * psa.r_d;
    @. g = FT(1e6) * g_lc;

    @. p_i = (g * p + a * r) / (a * k + g);
    @. psa.a_p = k * p_i;

    return nothing
);

product_limited_rate!(
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            apm::ApMethodC4VpmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack the variables from the cache
    f = cache.cache_incl_azi_2_1;
    qb = cache.cache_incl_azi_2_2;
    an = cache.cache_incl_azi_2_3;

    a = β * psa.v_pmax;
    d = psa.k_pep;
    p = air.s_aux.ps[2];
    r = β * psa.r_d;
    @. f = air.state.p_air / g_lc * FT(1e-6);

    qa = f;
    qc = a * p - r * (p + d);
    @. qb = f * r - p - d - a * f;
    @. an = lower_quadratic(qa, qb, qc);

    if g_lc[1] == 0 && g_lc[end] == 0
        @. psa.a_p = r;
    else
        @. psa.a_p = an + r;
    end;

    return nothing
);

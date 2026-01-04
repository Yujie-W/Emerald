"""

    product_limited_rate!(
                cache::SPACCache{FT},
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                g_lc::Vector{FT};
                β::FT = FT(1)) where {FT}

Update the product limited photosynthetic rate, given
- `cache` `SPACCache` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `p_i` Internal CO₂ partial pressure in `Pa`

"""
function product_limited_rate! end;

product_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate!(config, cache, ps.trait, ps.auxil, air, g_lc; β = β);

product_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT}, g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate_c3!(psa, config.METHODS.C3_AP_METHOD; β = β);

product_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT}, g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate_c4!(cache, pst, psa, config.METHODS.C4_AP_METHOD, air, g_lc; β = β);

product_limited_rate_c3!(
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC3Inf;
            β::FT = FT(1)) where {FT} = (@. psa.a_p = FT(Inf); return nothing);

product_limited_rate_c3!(
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC3Vcmax;
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.v_cmax / 2; return nothing);

product_limited_rate_c4!(
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VcmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack variables from the cache
    g = cache.cache_incl_azi_2_1;
    p_i = cache.cache_incl_azi_2_2;

    a = air.state.p_air;
    k = β * psa.k_pep_clm * pst.v_cmax25;
    p = air.auxil.ps[2];
    r = β * psa.r_d;
    @. g = FT(1e6) * g_lc;

    @. p_i = (g * p + a * r) / (a * k + g);
    @. psa.a_p = k * p_i;

    return nothing
);

product_limited_rate_c4!(
            cache::SPACCache{FT},
            ::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VpmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack the variables from the cache
    f = cache.cache_incl_azi_2_1;
    qb = cache.cache_incl_azi_2_2;
    an = cache.cache_incl_azi_2_3;

    a = β * psa.v_pmax;
    d = psa.k_pep;
    p = air.auxil.ps[2];
    r = β * psa.r_d;
    @. f = air.state.p_air / g_lc * FT(1e-6);

    qa = f;
    qc = a * p - r * (p + d);
    @. qb = f * r - p - d - a * f;
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_p[i] = r;
        else
            psa.a_p[i] = an[i] + r;
        end;
    end;

    return nothing
);

# Pressure mode for A-Ci fiting
product_limited_rate!(
            config::SPACConfig{FT},
            ps::LeafPhotosystem{FT},
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate!(config, ps.trait, ps.auxil, p_i; β = β);

product_limited_rate!(
            config::SPACConfig{FT},
            ::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate_c3!(psa, config.METHODS.C3_AP_METHOD; β = β);

product_limited_rate!(
            config::SPACConfig{FT},
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate_c4!(pst, psa, config.METHODS.C4_AP_METHOD, p_i; β = β);

product_limited_rate_c4!(
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VcmaxPi,
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.k_pep_clm * pst.v_cmax25 * p_i; return nothing);

product_limited_rate_c4!(
            ::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VpmaxPi,
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.v_pmax * p_i / (p_i + psa.k_pep); return nothing);

"""

    rubisco_limited_rate!(
                config::SPACConfig{FT},
                cache::SPACCache{FT},
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                g_lc::Vector{FT};
                β::FT = FT(1)) where {FT}

Update the RubisCO limited photosynthetic rate, given
- `cache` `SPACCache` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `p_i` Internal CO₂ partial pressure in `Pa`

"""
function rubisco_limited_rate! end;

rubisco_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(config, cache, ps.trait, ps.auxil, air, g_lc; β = β);

rubisco_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate_c3!(config, cache, psa, config.METHODS.C3_AC_METHOD, air, g_lc; β = β);

rubisco_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate_c4!(psa, config.METHODS.C4_AC_METHOD; β = β);

rubisco_limited_rate_c3!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            psa::LeafPhotosystemAuxil{FT},
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
    p = air.auxil.ps[2];
    r = β * psa.r_d;
    @. f = air.state.p_air / g_lc * FT(1e-6);

    qa = f;
    @. qb = f * r - p - d - a * f;
    @. qc = a * p - b - r * (p + d);
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_c[i] = r;
        else
            psa.a_c[i] = an[i] + r;
        end;
    end;

    return nothing
);

rubisco_limited_rate_c4!(
            psa::LeafPhotosystemAuxil{FT},
            ::AcMethodC4Vcmax;
            β::FT = FT(1)) where {FT} = (@. psa.a_c = β * psa.v_cmax; return nothing);

# Pressure mode for A-Ci fiting
rubisco_limited_rate!(
            config::SPACConfig{FT},
            ps::LeafPhotosystem{FT},
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate!(config, ps.trait, ps.auxil, p_i; β = β);

rubisco_limited_rate!(
            config::SPACConfig{FT},
            ::GeneralC3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate_c3!(psa, config.METHODS.C3_AC_METHOD, p_i; β = β);

rubisco_limited_rate!(
            config::SPACConfig{FT},
            ::GeneralC4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = rubisco_limited_rate_c4!(psa, config.METHODS.C4_AC_METHOD; β = β);

rubisco_limited_rate_c3!(
            psa::LeafPhotosystemAuxil{FT},
            ::AcMethodC3VcmaxPi,
            p_i::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_c = β * psa.v_cmax * (p_i - psa.γ_star) / (p_i + psa.k_m); return nothing);

"""

    light_limited_rate!(
                cache::SPACCache{FT},
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                g_lc::Vector{FT};
                β::FT = FT(1)) where {FT}

Update the electron transport limited photosynthetic rate, given
- `cache` `SPACCache` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function light_limited_rate! end;

light_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            ps::LeafPhotosystem{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = light_limited_rate!(config, cache, ps.trait, ps.state, ps.auxil, air, g_lc; β = β);

light_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = light_limited_rate!(config, cache, pss, psa, config.METHODS.C3_AJ_METHOD, air, g_lc; β = β);

light_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = light_limited_rate!(config, cache, pss, psa, config.METHODS.C4_AJ_METHOD, air, g_lc; β = β);

light_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ajm::AjMethodC3JmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if psa.j[1] == 0 && psa.j[end] == 0
        @. psa.a_j = 0;

        return nothing
    end;

    # unpack the cache variables
    b  = cache.cache_incl_azi_2_1;
    f  = cache.cache_incl_azi_2_2;
    qa = cache.cache_incl_azi_2_3;
    qb = cache.cache_incl_azi_2_4;
    qc = cache.cache_incl_azi_2_5;
    an = cache.cache_incl_azi_2_6;

    a = psa.j;
    c = pss.EFF_1;
    d = pss.EFF_2 * psa.γ_star;
    p = air.auxil.ps[2];
    r = β * psa.r_d;
    @. b = psa.j * psa.γ_star;
    @. f = air.state.p_air / g_lc * FT(1e-6);

    @. qa = c * f;
    @. qb = c * f * r - c * p - d - a * f;
    @. qc = a * p - b - r * (c * p + d);
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_j[i] = r;
        else
            psa.a_j[i] = an[i] + r;
        end;
    end;

    return nothing
);

light_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ajm::AjMethodC3VqmaxPi,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if psa.j_psi[1] == 0 && psa.j_psi[end] == 0
        @. psa.a_j = 0;

        return nothing
    end;

    eff_a = 1 - psa.η_l / psa.η_c;
    eff_b = 1 / psa.η_c;
    eff_1 = eff_a * pss.EFF_1 + 3 * eff_b;
    eff_2 = eff_a * pss.EFF_2 + 7 * eff_b;

    # unpack the cache variables
    b  = cache.cache_incl_azi_2_1;
    f  = cache.cache_incl_azi_2_2;
    qa = cache.cache_incl_azi_2_3;
    qb = cache.cache_incl_azi_2_4;
    qc = cache.cache_incl_azi_2_5;
    an = cache.cache_incl_azi_2_6;

    a  = psa.j_psi;
    c  = eff_1;
    d  = eff_2 * psa.γ_star;
    p  = air.auxil.ps[2];
    r  = β * psa.r_d;
    @. b = psa.j_psi * psa.γ_star;
    @. f = air.state.p_air / g_lc * FT(1e-6);

    @. qa = c * f;
    @. qb = c * f * r - c * p - d - a * f;
    @. qc = a * p - b - r * (c * p + d);
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_j[i] = r;
        else
            psa.a_j[i] = an[i] + r;
        end;
    end;

    return nothing
);

light_limited_rate!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ajm::AjMethodC4JPSII,
            air::AirLayer{FT},
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (@. psa.a_j = psa.j * psa.e2c; return nothing);

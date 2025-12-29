"""

    photosystem_electron_transport!(
                config::SPACConfig{FT},
                cache::SPACCache{FT},
                ps::CanopyLayerPhotosystem{FT},
                ppar::Vector{FT},
                p_i::Union{FT, Vector{FT}};
                β::FT = FT(1)) where {FT}
    photosystem_electron_transport!(
                config::SPACConfig{FT},
                ps::LeafPhotosystem{FT},
                ppar::FT,
                p_i::FT;
                β::FT = FT(1)) where {FT}

Update the electron transport rates, given
- `cache` `SPACCache` type struct
- `ps` `CanopyLayerPhotosystem` or `LeafPhotosystem` type struct
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_electron_transport! end;

# For Leaf
photosystem_electron_transport!(
            config::SPACConfig{FT},
            ps::LeafPhotosystem{FT},
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(config.METHODS.AJM, config.METHODS.COLIMIT_J, ps.trait, ps.state, ps.auxil, ppar, p_i; β = β);

photosystem_electron_transport!(
            ::AjMethodC3JmaxPi,
            colimj::UnionColimit{FT},
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (
    psa.e2c   = (p_i == Inf) ? (1 / pss.EFF_1) : (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    psa.j     = colimited_rate(psa.j_pot, β * psa.j_max, colimj);

    return nothing
);

photosystem_electron_transport!(
            ::AjMethodC3VqmaxPi,
            colimj::SerialColimit,
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (
    psa.e2c   = (p_i == Inf) ? (1 / pss.EFF_1) : (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    psa.j_psi = colimited_rate(β * psa.v_qmax, ppar * (1 - psa.f_psii) * psa.ϕ_psi_max, colimj);
    psa.η     = (p_i == Inf) ? (1 - psa.η_l / psa.η_c + 3 / pss.EFF_1 / psa.η_c) : (1 - psa.η_l / psa.η_c + (3 * p_i + 7 * psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star) / psa.η_c);
    psa.j_pot = psa.j_psi / psa.η;
    psa.j     = psa.j_pot;

    return nothing
);

photosystem_electron_transport!(
            ::AjMethodC4JPSII,
            ::UnionColimit{FT},
            pst::GeneralC4Trait{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (
    psa.e2c   = 1 / 6;
    psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    psa.j     = psa.j_pot;

    return nothing
);

# For CanopyLayer
photosystem_electron_transport!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            ps::CanopyLayerPhotosystem{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(config, cache, ps.state, ps.auxil, ppar, p_i; β = β);

photosystem_electron_transport!(
            config::SPACConfig{FT},
            cache::SPACCache{FT},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport_c3!(config.METHODS.C3_AJ_METHOD, config.METHODS.COLIMIT_J, cache, pss, psa, ppar, p_i; β = β);

photosystem_electron_transport!(
            config::SPACConfig{FT},
            ::SPACCache{FT},
            pss::C4State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport_c4!(config.METHODS.C4_AJ_METHOD, psa, ppar);

photosystem_electron_transport_c3!(
            ::AjMethodC3JmaxPi,
            colimj::UnionColimit{FT},
            cache::SPACCache{FT},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    @. psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    @. psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    colimited_rate!(β * psa.j_max, psa.j_pot, psa.j, colimj);
    for i in eachindex(psa.e2c)
        isnan(psa.e2c[i]) ? (psa.e2c[i] = 1 / pss.EFF_1) : nothing;
    end;

    return nothing
);

photosystem_electron_transport_c3!(
            ::AjMethodC3VqmaxPi,
            colimj::SerialColimit,
            cache::SPACCache{FT},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    _j = cache.cache_incl_azi_2_1;
    @. _j = ppar * (1 - psa.f_psii) * psa.ϕ_psi_max;
    colimited_rate!(β * psa.v_qmax, _j, psa.j_psi, colimj);

    @. psa.η = 1 - psa.η_l / psa.η_c + (3 * p_i + 7 * psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star) / psa.η_c;
    for i in eachindex(psa.η)
        isnan(psa.η[i]) ? (psa.η[i] = 1 - psa.η_l / psa.η_c + 3 / pss.EFF_1 / psa.η_c) : nothing;
    end;
    @. psa.j_pot = psa.j_psi / psa.η;
    @. psa.j = psa.j_pot;
    @. psa.e2c = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    for i in eachindex(psa.e2c)
        isnan(psa.e2c[i]) ? (psa.e2c[i] = 1 / pss.EFF_1) : nothing;
    end;

    return nothing
);

photosystem_electron_transport_c4!(::AjMethodC4JPSII, psa::CanopyLayerPhotosystemAuxil{FT}, ppar::Vector{FT}) where {FT} = (
    @. psa.e2c   = 1 / 6;
    @. psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    @. psa.j     = psa.j_pot;

    return nothing
);

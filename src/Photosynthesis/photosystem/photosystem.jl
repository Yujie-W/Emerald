"""

    photosynthesis!(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; use_glc::Bool = false) where {FT}
    photosynthesis!(config::SPACConfig{FT}, cache::SPACCache{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::Vector{FT}, ppar::Vector{FT}, t::FT) where {FT}

Update the photosynthesis rate, given
- `config` `SPACConfig` struct
- `cache` `SPACCache` struct
- `leaf` `Leaf` struct
- `air` `AirLayer` struct
- `use_glc` Whether to use leaf conductance to CO₂ to compute photosynthesis rate (otherwise use intercellular CO₂ partial pressure)
- `ps` `LeafPhotosystem` struct
- `p_i` Vector of intercellular CO₂ partial pressures
- `ppar` Vector of photosynthetically active radiation
- `t` Leaf temperature in Kelvin

"""
function photosynthesis! end;

photosynthesis!(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; use_glc::Bool = false) where {FT} = (
    if use_glc
        photosynthesis_glc!(config, cache, leaf, air);
    else
        photosynthesis!(config, cache, leaf.photosystem, air, leaf.flux.auxil.p_CO₂_i, leaf.flux.auxil.ppar, leaf.energy.auxil.t);
    end;

    return nothing
);

photosynthesis!(config::SPACConfig{FT}, cache::SPACCache{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::Vector{FT}, ppar::Vector{FT}, t::FT) where {FT} = (
    photosystem_temperature_dependence!(config, ps, air, t);
    photosystem_electron_transport!(config, cache, ps, ppar, p_i);
    rubisco_limited_rate!(config, ps, p_i);
    light_limited_rate!(ps);
    product_limited_rate!(config, ps, p_i);
    colimit_photosynthesis!(config, ps);
    photosystem_coefficients!(config, cache, ps, ppar);

    return nothing
);

photosynthesis_glc!(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    ps = leaf.photosystem;
    g_lc = leaf.flux.auxil.g_CO₂;
    ppar = leaf.flux.auxil.ppar;
    t = leaf.energy.auxil.t;
    p_i = leaf.flux.auxil.p_CO₂_i;

    photosystem_temperature_dependence!(config, ps, air, t);
    photosystem_electron_transport!(config, cache, ps, ppar, p_i);
    rubisco_limited_rate!(config, cache, ps, air, g_lc);
    light_limited_rate!(config, cache, ps, air, g_lc);
    product_limited_rate!(config, cache, ps, air, g_lc);
    colimit_photosynthesis!(config, ps);
    p_i = @. air.auxil.ps[2] - ps.auxil.a_n / g_lc * air.state.p_air .* FT(1e-6);
    photosystem_electron_transport!(config, cache, ps, ppar, p_i);
    photosystem_coefficients!(config, cache, ps, ppar);

    return nothing
);

"""

    photosynthesis!(config::SPACConfig{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::FT, ppar::FT, t::FT) where {FT}

Update the photosynthesis rate, given
- `config` `SPACConfig` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `ppar` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `t` Leaf temperature in `K`

"""
function photosynthesis! end;

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

"""

    aci_an(config::SPACConfig{FT}, cache::SPACCache{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::Number, ppar::Number, t::Number) where {FT}

Compute the net photosynthetic rate, given
- `config` `SPACConfig` struct
- `cache` `SPACCache` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `ppar` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `t` Leaf temperature in `K`

"""
function aci_an(config::SPACConfig{FT}, cache::SPACCache{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::FT, ppar::FT, t::FT) where {FT}
    photosynthesis!(config, cache, ps, air, [p_i,], [ppar,], t);

    return ps.auxil.a_n[1]
end;


"""

    aci_curve(config::SPACConfig{FT},
              cache::SPACCache{FT},
              ps::LeafPhotosystem{FT},
              air::AirLayer{FT},
              pis::Vector,
              ppars::Vector,
              ts::Vector) where {FT}
    aci_curve(config::SPACConfig{FT},
              cache::SPACCache{FT},
              ps::LeafPhotosystem{FT},
              air::AirLayer{FT},
              df::DataFrame) where {FT}

Compute the net photosynthetic rates, given
- `config` `SPACConfig` struct
- `cache` `SPACCache` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `pis` Internal CO₂ partial pressure in `Pa`
- `ppars` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `ts` Leaf temperature in `K`
- `df` DataFrame with columns `P_I`, `PPAR`, and `T_LEAF`

"""
function aci_curve end;

aci_curve(config::SPACConfig{FT},
          cache::SPACCache{FT},
          ps::LeafPhotosystem{FT},
          air::AirLayer{FT},
          pis::Vector,
          ppars::Vector,
          ts::Vector) where {FT} = aci_an.((config,), (cache,), (ps,), (air,), pis, ppars, ts);

aci_curve(config::SPACConfig{FT},
          cache::SPACCache{FT},
          ps::LeafPhotosystem{FT},
          air::AirLayer{FT},
          df::DataFrame) where {FT} = aci_curve(config, cache, ps, air, df.P_I, df.PPAR, df.T_LEAF);

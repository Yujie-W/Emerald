#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain net photosynthetic rates of given data
#
#######################################################################################################################################################################################################
"""

    aci_curve(config::SPACConfig{FT},
              ps::LeafPhotosystem{FT},
              air::AirLayer{FT},
              pis::Vector,
              ppars::Vector,
              ts::Vector) where {FT}
    aci_curve(config::SPACConfig{FT},
              ps::LeafPhotosystem{FT},
              air::AirLayer{FT},
              df::DataFrame) where {FT}

Compute the net photosynthetic rates, given
- `config` `SPACConfig` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `pis` Internal CO₂ partial pressure in `Pa`
- `ppars` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `ts` Leaf temperature in `K`
- `df` DataFrame with columns `P_I`, `PPAR`, and `T_LEAF`

"""
function aci_curve end;

aci_curve(config::SPACConfig{FT},
          ps::LeafPhotosystem{FT},
          air::AirLayer{FT},
          pis::Vector,
          ppars::Vector,
          ts::Vector) where {FT} = aci_an.((config,), (ps,), (air,), pis, ppars, ts);

aci_curve(config::SPACConfig{FT},
          ps::LeafPhotosystem{FT},
          air::AirLayer{FT},
          df::DataFrame) where {FT} = aci_curve(config, ps, air, df.P_I, df.PPAR, df.T_LEAF);

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain net photosynthetic rate
#
#######################################################################################################################################################################################################
"""

    aci_an(config::SPACConfiguration{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::Number, ppar::Number, t::Number) where {FT}

Compute the net photosynthetic rate, given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `ppar` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `t` Leaf temperature in `K`

"""
function aci_an(config::SPACConfiguration{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::Number, ppar::Number, t::Number) where {FT}
    photosynthesis!(config, ps, air, p_i, ppar, t);

    return ps.auxil.a_n
end;

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add photosynthesis! function to compute photosynthesis rate from Pi, PPAR, and Tleaf
#
#######################################################################################################################################################################################################
"""

    photosynthesis!(ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::FT, ppar::FT, t::FT) where {FT}

Update the photosynthesis rate, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `ppar` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `t` Leaf temperature in `K`

"""
function photosynthesis! end;

photosynthesis!(config::SPACConfiguration{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::FT, ppar::FT, t::FT) where {FT} = (
    photosystem_temperature_dependence!(config, ps, air, t);
    photosystem_electron_transport!(ps, ppar, p_i);
    rubisco_limited_rate!(ps, p_i);
    light_limited_rate!(ps);
    product_limited_rate!(ps, p_i);
    colimit_photosynthesis!(ps);
    photosystem_coefficients!(config, ps, ppar);

    return nothing
);

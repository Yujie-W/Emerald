# This file contains function to compute photosynthetic rates only (to use with optimality model)

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method to compute photosynthetic rates only
#     2024-Jul-24: add new method with TD to speed up (need to call shaded part first to update TD)
#
#######################################################################################################################################################################################################
"""

    photosynthesis_only!(psm::Union{C3Cyto{FT}, C3VJP{FT}, C4VJP{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT) where {FT}

Updates leaf photosynthetic rates based on leaf diffusive conductance (for StomataModels.jl temporary use), given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` type photosynthesis system
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `ppar` PPAR used for photosynthesis
- `t` Leaf temperature in `[K]`

"""
function photosynthesis_only! end;

photosynthesis_only!(psm::CanopyLayerPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT) where {FT} = (
    #=
    photosystem_temperature_dependence!(psm, air, t);
    photosystem_electron_transport!(psm, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(psm, air, g_lc; β = FT(1));
    light_limited_rate!(psm, air, g_lc; β = FT(1));
    product_limited_rate!(psm, air, g_lc; β = FT(1));
    colimit_photosynthesis!(psm; β = FT(1));
    =#

    return FT(1)
);

photosynthesis_only!(psm::CanopyLayerPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT, ppar::FT) where {FT} = (
    #=
    photosystem_electron_transport!(psm, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(psm, air, g_lc; β = FT(1));
    light_limited_rate!(psm, air, g_lc; β = FT(1));
    product_limited_rate!(psm, air, g_lc; β = FT(1));
    colimit_photosynthesis!(psm; β = FT(1));
    =#

    return FT(1)
);

photosynthesis_only!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT) where {FT} = (
    photosystem_temperature_dependence!(psm, air, t);
    photosystem_electron_transport!(psm, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(psm, air, g_lc; β = FT(1));
    light_limited_rate!(psm, air, g_lc; β = FT(1));
    product_limited_rate!(psm, air, g_lc; β = FT(1));
    colimit_photosynthesis!(psm; β = FT(1));

    return psm.auxil.a_n
);

photosynthesis_only!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::FT, ppar::FT) where {FT} = (
    photosystem_electron_transport!(psm, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(psm, air, g_lc; β = FT(1));
    light_limited_rate!(psm, air, g_lc; β = FT(1));
    product_limited_rate!(psm, air, g_lc; β = FT(1));
    colimit_photosynthesis!(psm; β = FT(1));

    return psm.auxil.a_n
);

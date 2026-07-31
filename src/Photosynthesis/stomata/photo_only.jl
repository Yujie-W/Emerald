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

    photosynthesis_only!(config::SPACConfig{FT}, psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::Vector{FT}, ppar::Vector{FT}) where {FT}

Updates leaf photosynthetic rates based on leaf diffusive conductance (for StomataModels.jl temporary use), given
- `config` `SPACConfig` type structure
- `psm` `LeafPhotosystem` type structure
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `ppar` PPAR used for photosynthesis

"""
function photosynthesis_only! end;

photosynthesis_only!(config::SPACConfig{FT}, psm::LeafPhotosystem{FT}, air::AirLayer{FT}, g_lc::Vector{FT}, ppar::Vector{FT}) where {FT} = (
    (; NEW_C4_STOMATAL_METHODS) = config.FEATURES;

    photosystem_electron_transport!(config.METHODS.PS_METHODS, psm, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(config.METHODS.PS_METHODS, psm, air.state.p_air, air.auxil.ps[2], g_lc; β = FT(1));
    light_limited_rate!(config.METHODS.PS_METHODS, psm, air.state.p_air, air.auxil.ps[2], g_lc; β = FT(1));
    product_limited_rate!(config.METHODS.PS_METHODS, psm, air.state.p_air, air.auxil.ps[2], g_lc; β = FT(1));
    colimit_photosynthesis!(config.METHODS.PS_METHODS, psm; β = FT(1));

    return NEW_C4_STOMATAL_METHODS ? min.(psm.auxil.a_p, psm.auxil.a_j) .- psm.auxil.r_d : psm.auxil.a_n
);

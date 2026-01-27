"""

    leaf_level_spac_cache(config::SPACConfig{FT}) where {FT}

Create and return a `SPACCache` struct for leaf-level simulations, given
- `config` `SPACConfig` struct

"""
function leaf_level_spac_cache(config::SPACConfig{FT}) where {FT}
    return SPACCache{FT}(
                config.DIMENSIONS.DIM_AZI,
                config.DIMENSIONS.DIM_INCL,
                1,
                0,
                length(config.CONSTANTS.SPECTRA.Λ_SIF),
                length(config.CONSTANTS.SPECTRA.Λ_SIFE),
                length(config.CONSTANTS.SPECTRA.Λ))
end;

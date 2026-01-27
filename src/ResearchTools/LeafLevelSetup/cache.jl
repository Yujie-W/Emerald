"""

    leaf_level_spac_cache(FT)
    leaf_level_spac_cache(config::SPACConfig{FT}) where {FT}

Create and return a `SPACCache` struct for leaf-level simulations, given
- `FT` Floating-point type (e.g., `Float32`, `Float64`)
- `config` `SPACConfig` struct

"""
function leaf_level_spac_cache end;

leaf_level_spac_cache(FT) = leaf_level_spac_cache(leaf_level_config(FT));

leaf_level_spac_cache(config::SPACConfig{FT}) where {FT} = (
    @assert config.DIMENSIONS.DIM_PPAR_BINS == 0 "For leaf-level simulations, DIM_PPAR_BINS must be set to 0!";

    return SPACCache(config, 1)
);

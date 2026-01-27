"""

    leaf_level_config(FT = Float64)

Create and return a `SPACConfig` struct for leaf-level simulations, given
- `FT` Floating-point type (e.g., `Float32`, `Float64`)

"""
function leaf_level_config(FT = Float64)
    config = SPACConfig{FT}();
    config.DIMENSIONS.DIM_PPAR_BINS = 0;

    return config
end;

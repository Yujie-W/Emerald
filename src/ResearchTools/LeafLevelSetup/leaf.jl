"""

    leaf_level_leaf(FT, c3c4::String)
    leaf_level_leaf(config::SPACConfig{FT}, c3c4::String) where {FT}

Create and return a `Leaf` struct for leaf-level simulations, given
- `FT` Floating-point type (e.g., `Float32`, `Float64`)
- `c3c4` String indicating whether to create a C3 or C4 leaf
- `config` `SPACConfig` struct

"""
function leaf_level_leaf end;

leaf_level_leaf(FT, c3c4::String) = leaf_level_leaf(leaf_level_config(FT), c3c4);

leaf_level_leaf(config::SPACConfig{FT}, c3c4::String) where {FT} = (
    @assert c3c4 in ["C3", "C4"] "The model string should be either C3 or C4!";

    return Leaf(config, c3c4)
);

"""

    ∂A∂E_∂Θ∂E(config::SPACConfig{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}
    ∂A∂E_∂Θ∂E(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}

Return the ∂A∂E and ∂Θ∂E, given
- `config` `SPACConfig` struct
- `leaf` `Leaf` struct
- `air` `AirLayer` struct
- `cache` `SPACCache` struct (optional, can be generated from `config`)

"""
function ∂A∂E_∂Θ∂E end;

∂A∂E_∂Θ∂E(config::SPACConfig{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = ∂A∂E_∂Θ∂E(config, leaf_level_spac_cache(config), leaf, air);

∂A∂E_∂Θ∂E(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    leaf.flux.auxil.g_CO₂[1] = 1 / (1 / leaf.flux.auxil.g_CO₂_b + 1.6 / leaf.flux.state.g_H₂O_s[1] + 1 / leaf.flux.auxil.g_m[1]);

    # Calculate total conductance (g) and leaf-to-air VPD (d), and then update leaf xylem flux
    g_h = 1 / (1 / leaf.flux.state.g_H₂O_s[1] + 1 / (1.35 * leaf.flux.auxil.g_CO₂_b[1]));
    d = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.state.p_leaf * 1e6) - air.auxil.ps[3];
    f = g_h * d / air.state.p_air * leaf.xylem.trait.area;
    set_flow_profile!(leaf.xylem, f - leaf.capacitor.auxil.flow);

    # compute the dg/dt
    leaf_pressure_profile!(config, leaf, cache, leaf.xylem.auxil.pressure[1]);
    leaf_photosynthesis!(config, cache, leaf, air);
    ∂g∂t!(config, cache, leaf, air);

    return leaf.flux.auxil.∂A∂E[1], leaf.flux.auxil.∂Θ∂E[1]
);

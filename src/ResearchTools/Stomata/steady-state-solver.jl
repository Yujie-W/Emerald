"""

    steady_state_gs!(config::SPACConfig{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; timer::Number = 50000) where {FT}
    steady_state_gs!(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; timer::Number = 50000) where {FT}

Compute the steady-state stomatal conductance for a leaf, given
- `config` `SPACConfig` struct
- `leaf` `Leaf` struct
- `air` `AirLayer` struct
- `timer` Maximum time allowed for reaching steady state `[s]`
- `cache` `SPACCache` struct (optional, can be generated from `config`)

"""
function steady_state_gs! end;

steady_state_gs!(config::SPACConfig{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; timer::Number = 50000) where {FT} =
    steady_state_gs!(config, leaf_level_spac_cache(config), leaf, air; timer = timer);

steady_state_gs!(config::SPACConfig{FT}, cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; timer::Number = 50000) where {FT} = (
    @assert config.DIMENSIONS.DIM_PPAR_BINS == 0 "steady_state_gs! only supports leaf-level simulations (DIM_PPAR_BINS == 0)";

    δt_remain::FT = timer;
    while δt_remain > 0
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

        if abs(leaf.flux.auxil.∂g∂t[1]) <= 1e-7
            break
        end;

        # adjust time step based on capacitance buffer and stomatal conductance change rate
        δt = dynamic_timer(leaf, δt_remain);
        stomatal_conductance!(leaf, δt)
        leaf_water_budget!(leaf, leaf.xylem.auxil, δt);
        substep_aux!(leaf, false);
        δt_remain -= δt;
    end;

    return nothing
);


"""

    dynamic_timer(leaf::Leaf{FT}, δt::FT) where {FT}

Adjust the time step based on leaf capacitance buffer and stomatal conductance change rate, given
- `leaf` `Leaf` struct
- `δt` Current time step `[s]`

"""
function dynamic_timer(leaf::Leaf{FT}, δt::FT) where {FT}
    new_δt::FT = min(10, δt);

    if leaf.capacitor.auxil.flow > 0
        new_v = capacitance_volume(leaf.capacitor.trait.pv, leaf.capacitor.auxil.p - FT(0.1), leaf.energy.auxil.t) * leaf.capacitor.trait.v_max * leaf.xylem.trait.area;
        new_δt = min((leaf.capacitor.state.v_storage - new_v) / leaf.capacitor.auxil.flow, new_δt);
    elseif leaf.capacitor.auxil.flow < 0
        new_v = capacitance_volume(leaf.capacitor.trait.pv, leaf.capacitor.auxil.p + FT(0.1), leaf.energy.auxil.t) * leaf.capacitor.trait.v_max * leaf.xylem.trait.area;
        new_δt = min((leaf.capacitor.state.v_storage - new_v) / leaf.capacitor.auxil.flow, new_δt);
    end;
    if isnan(new_δt)
        @error "NaN or very small δt detected when adjusting δt based on leaf capacitance buffer" leaf.capacitor.auxil.p leaf.capacitor.auxil.flow;
        return error("NaN detected in dynamic_timer")
    end;

    # make sure each leaf stomatal conductances do not change more than 0.001 mol m⁻² s⁻¹
    for ∂g∂t in leaf.flux.auxil.∂g∂t
        new_δt = min(FT(0.001) / abs(∂g∂t), new_δt);
        if isnan(new_δt)
            @error "NaN or very small δt detected when adjusting δt based on leaf stomatal conductance dYdt" ∂g∂t;
            return error("NaN detected in dynamic_timer")
        end;
    end;

    return new_δt
end;

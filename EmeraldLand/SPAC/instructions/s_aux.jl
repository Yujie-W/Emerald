#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-26: add s_aux! method for AirLayerState-dependent variables
#     2024-Feb-26: add s_aux! method for AirLayer
#     2024-Feb-27: add s_aux! method for SoilLayerState-dependent variables
#     2024-Feb-27: add s_aux! method for SoilLayer
#     2024-Feb-27: add s_aux! method for Root
#     2024-Feb-27: add s_aux! method for JunctionCapacitor
#     2024-Feb-27: add s_aux! method for Stem
#     2024-Feb-27: add s_aux! method for Leaf
#     2024-Feb-27: add s_aux! method for BulkSPAC
#
#######################################################################################################################################################################################################
s_aux!(t_aux::AirLayerTDAuxil{FT}, state::AirLayerState{FT}, s_aux::AirLayerSDAuxil{FT}) where {FT} = (
    s_aux.t = state.Σe / heat_capacitance(state);
    for i in 1:5
        s_aux.ps[i] = (state.ns[i] * GAS_R(FT) * s_aux.t) / t_aux.δz;
    end;
    s_aux.f_CO₂ = s_aux.ps[2] / state.p_air * 1e6;

    return nothing
);

s_aux!(air::AirLayer{FT}) where {FT} = s_aux!(air.t_aux, air.state, air.s_aux);

s_aux!(trait::SoilLayerTrait{FT}, t_aux::SoilLayerTDAuxil{FT}, state::SoilLayerState{FT}, s_aux::SoilLayerSDAuxil{FT}) where {FT} = (
    s_aux.cp = heat_capacitance(trait, t_aux, state);
    s_aux.t = state.Σe / s_aux.cp;

    # update the conductance, potential, diffusivity, and thermal conductivity (0.5 for tortuosity factor)
    s_aux.k = relative_soil_k(trait.vc, state.θ) * trait.vc.K_MAX * relative_viscosity(s_aux.t) / t_aux.δz;
    s_aux.ψ = soil_ψ_25(trait.vc, state.θ; oversaturation = true) * relative_surface_tension(s_aux.t);
    s_aux.kd = 0.5 * max(0, trait.vc.Θ_SAT - state.θ) / t_aux.δz;
    s_aux.kv = 0.5 * trait.vc.Θ_SAT / max(FT(0.01), trait.vc.Θ_SAT - state.θ) / t_aux.δz;
    s_aux.λ_soil_water = (trait.λ_soil + state.θ * Λ_THERMAL_H₂O(FT)) / t_aux.δz;

    return nothing
);

s_aux!(soil::SoilLayer{FT}) where {FT} = s_aux!(soil.trait, soil.t_aux, soil.state, soil.s_aux);

s_aux!(root::Root{FT}) where {FT} = (
    # update root cp and temperature
    root.energy.s_aux.cp = heat_capacitance(root);
    root.energy.s_aux.t = root.energy.state.Σe / root.energy.s_aux.cp;

    return nothing
);

s_aux!(junc::JunctionCapacitor{FT}) where {FT} = (
    # update junction cp and temperature
    junc.s_aux.cp = heat_capacitance(junc);
    junc.s_aux.t = junc.state.Σe / junc.s_aux.cp;

    # update the junction buffer pressure
    junc.s_aux.pressure = capacitance_pressure(junc.trait.pv, junc.state.v_storage / junc.trait.v_max, junc.s_aux.t);

    return nothing
);

s_aux!(stem::Stem{FT}) where {FT} = (
    # update stem cp and temperature
    stem.energy.s_aux.cp = heat_capacitance(stem);
    stem.energy.s_aux.t = stem.energy.state.Σe / stem.energy.s_aux.cp;

    return nothing
);

s_aux!(leaf::Leaf{FT}) where {FT} = (
    s_aux!(leaf.capacitor.state, leaf.xylem.trait, leaf.bio.trait, leaf.energy.state, leaf.energy.s_aux);

    return nothing
);

s_aux!(plant::Plant{FT}) where {FT} = (
    # update the s_aux for each of the field in the plant
    for root in plant.roots
        s_aux!(root);
    end;
    s_aux!(plant.junction);
    s_aux!(plant.trunk);
    for stem in plant.branches
        s_aux!(stem);
    end;
    for leaf in plant.leaves
        s_aux!(leaf);
    end;

    return nothing
);

s_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    # update the s_aux for each of the field in the bulk spac system (the order might matter, will need to check)
    # the soil auxilary variables
    for soil in spac.soils
        s_aux!(soil);
    end;

    # the air auxilary variables
    for air in spac.airs
        s_aux!(air);
    end;

    # the plant auxilary variables
    s_aux!(spac.plant);

    # the canopy auxilary variables
    s_aux!(config, spac.canopy);

    return nothing
);

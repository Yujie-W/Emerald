#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-26: add s_aux! method for AirLayerState-dependent variables
#     2024-Feb-26: add s_aux! method for AirLayer
#     2024-Feb-26: add method for LeafEnergyState-dependent auxiliary variables
#     2024-Feb-27: add s_aux! method for SoilLayerState-dependent variables
#     2024-Feb-27: add s_aux! method for SoilLayer
#     2024-Feb-27: add s_aux! method for Root
#     2024-Feb-27: add s_aux! method for JunctionCapacitor
#     2024-Feb-27: add s_aux! method for Stem
#     2024-Feb-27: add s_aux! method for Leaf
#     2024-Feb-27: add s_aux! method for BulkSPAC
#     2024-Jul-24: add leaf shedded flag
#     2024-Jul-30: compute OCS fraction in the air layer
#     2024-Nov-05: remove leaf shedded flag
#     2025-Jun-05: make soil total energy relative to triple temperature for phase change purposes
#
#######################################################################################################################################################################################################
"""

    s_aux!(spac::BulkSPAC{FT}) where {FT}

Update the prognostic state variable (such as total energy, water storage, etc.) dependent auxiliary variables for the SPAC system, given
- `spac` SPAC

"""
function s_aux! end;

s_aux!(spac::BulkSPAC{FT}) where {FT} = (
    # update the s_aux for each of the field in the bulk spac system (the order might matter, will need to check)
    # the soil auxiliary variables
    for soil in spac.soils
        s_aux!(soil);
    end;

    # the air auxiliary variables
    for air in spac.airs
        s_aux!(air);
    end;

    # the plant auxiliary variables
    s_aux!(spac.plant);

    return nothing
);

s_aux!(soil::SoilLayer{FT}) where {FT} = (
    soil.s_aux.cp = heat_capacitance(soil);
    soil.s_aux.t = soil.state.Σe / soil.s_aux.cp + T₀(FT);

    # update the conductance, potential, diffusivity, and thermal conductivity (0.5 for tortuosity factor)
    soil.s_aux.k = relative_soil_k(soil.trait.vc, soil.state.θ) * soil.trait.vc.K_MAX * relative_viscosity(soil.s_aux.t) / soil.t_aux.δz;
    soil.s_aux.ψ = soil_ψ_25(soil.trait.vc, soil.state.θ; oversaturation = true) * relative_surface_tension(soil.s_aux.t);
    soil.s_aux.kd = 0.5 * max(0, soil.trait.vc.Θ_SAT - soil.state.θ) / soil.t_aux.δz;
    soil.s_aux.kv = 0.5 * soil.trait.vc.Θ_SAT / max(FT(0.01), soil.trait.vc.Θ_SAT - soil.state.θ) / soil.t_aux.δz;
    soil.s_aux.λ_soil_water = (soil.trait.λ_soil + soil.state.θ * Λ_THERMAL_H₂O(FT)) / soil.t_aux.δz;

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

s_aux!(leaf::CanopyLayer{FT}) where {FT} = (
    if leaf.xylem.trait.area > 0
        leaf.energy.s_aux.cp = heat_capacitance(leaf);
        leaf.energy.s_aux.t = leaf.energy.state.Σe / leaf.energy.s_aux.cp;
    end;

    return nothing
);

s_aux!(air::AirLayer{FT}) where {FT} = (
    air.s_aux.t = air.state.Σe / heat_capacitance(air);
    for i in 1:6
        air.s_aux.ps[i] = (air.state.ns[i] * GAS_R(FT) * air.s_aux.t) / air.t_aux.δz;
    end;
    air.s_aux.f_CO₂ = air.s_aux.ps[2] / air.state.p_air * 1e6;
    air.s_aux.f_OCS = air.s_aux.ps[6] / air.state.p_air * 1e9;

    return nothing
);

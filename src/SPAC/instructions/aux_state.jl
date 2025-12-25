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
#     2025-Jun-05: account for ice volume in the calculation of vapor and gas diffusion related auxiliary variables
#     2025-Jun-07: set a maximum effective water volume to each layer as the soil is allowed to be oversaturated
#
#######################################################################################################################################################################################################
"""

    s_aux!(spac::BulkSPAC{FT}) where {FT}

Update the prognostic state variable (such as total energy, water storage, etc.) dependent auxiliary variables for the SPAC system, given
- `spac` SPAC

"""
function s_aux! end;

s_aux!(spac::BulkSPAC{FT}) where {FT} = (
    # update the state-dependent auxil for each of the field in the bulk spac system (the order might matter, will need to check)
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
    soil.auxil.cp = heat_capacitance(soil);
    soil.auxil.t = soil.state.Σe / soil.auxil.cp + T₀(FT);

    # update the conductance, potential, diffusivity, and thermal conductivity (0.5 for tortuosity factor)
    # TODO: add Λ_THERMAL_H₂O for ice
    soil.auxil.k = relative_soil_k(soil.trait.vc, soil.state.θ) * soil.trait.vc.K_MAX * relative_viscosity(soil.auxil.t) / soil.auxil.δz;
    soil.auxil.ψ = soil_ψ_25(soil.trait.vc, soil.state.θ; oversaturation = true) * relative_surface_tension(soil.auxil.t);
    soil.auxil.kd = 0.5 * max(0, soil.trait.vc.Θ_SAT - soil.state.θ - soil.state.θ_ice) / soil.auxil.δz;
    soil.auxil.kv = 0.5 * soil.trait.vc.Θ_SAT / max(FT(0.01), soil.trait.vc.Θ_SAT - soil.state.θ - soil.state.θ_ice) / soil.auxil.δz;
    soil.auxil.λ_soil_water = (soil.trait.λ_soil + max(soil.trait.vc.Θ_SAT, soil.state.θ + soil.state.θ_ice) * Λ_THERMAL_H₂O(FT)) / soil.auxil.δz;

    return nothing
);

s_aux!(plant::Plant{FT}) where {FT} = (
    # update the state-dependent auxil for each of the field in the plant
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
    root.energy.auxil.cp = heat_capacitance(root);
    root.energy.auxil.t = root.energy.state.Σe / root.energy.auxil.cp;

    return nothing
);

s_aux!(junc::JunctionCapacitor{FT}) where {FT} = (
    # update junction cp and temperature
    junc.auxil.cp = heat_capacitance(junc);
    junc.auxil.t = junc.state.Σe / junc.auxil.cp;

    # update the junction buffer pressure
    junc.auxil.pressure = capacitance_pressure(junc.trait.pv, junc.state.v_storage / junc.trait.v_max, junc.auxil.t);

    return nothing
);

s_aux!(stem::Stem{FT}) where {FT} = (
    # update stem cp and temperature
    stem.energy.auxil.cp = heat_capacitance(stem);
    stem.energy.auxil.t = stem.energy.state.Σe / stem.energy.auxil.cp;

    return nothing
);

s_aux!(leaf::CanopyLayer{FT}) where {FT} = (
    if leaf.xylem.trait.area > 0
        leaf.energy.auxil.cp = heat_capacitance(leaf);
        leaf.energy.auxil.t = leaf.energy.state.Σe / leaf.energy.auxil.cp;
    end;

    return nothing
);

s_aux!(air::AirLayer{FT}) where {FT} = (
    air.auxil.t = air.state.Σe / heat_capacitance(air);
    for i in 1:6
        air.auxil.ps[i] = (air.state.ns[i] * GAS_R(FT) * air.auxil.t) / air.auxil.δz;
    end;
    air.auxil.f_CO₂ = air.auxil.ps[2] / air.state.p_air * 1e6;
    air.auxil.f_OCS = air.auxil.ps[6] / air.state.p_air * 1e9;

    return nothing
);

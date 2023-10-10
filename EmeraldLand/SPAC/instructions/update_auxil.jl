# This file contains functions to update the auxiliary variables at different time steps
#     - Big time step (such as the integrated values: GPP, ET, etc.)
#     - Sub time step within a big time step (such as the cache variables which can be recomputed: t, cp, and p)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_substep_auxils! function
#     2023-Sep-30: add buffer pressure and flow calculations
#     2023-Oct-06: add soil auxiliary variable calculations
#     2023-Oct-07: add soil bulk auxiliary variable calculations
#
#######################################################################################################################################################################################################
"""

    update_substep_auxils!(spac::MultiLayerSPAC{FT}) where {FT}

Update the auxiliary variables at sub time step within a big time step, given
- `spac` `MultiLayerSPAC` SPAC

"""
function update_substep_auxils! end;

update_substep_auxils!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; SOILS, SOIL_BULK, ROOTS, JUNCTION, TRUNK, BRANCHES, LEAVES) = spac;

    # update the soil layer auxiliary variables
    for soil in SOILS
        update_substep_auxils!(soil);
    end;

    # update soil bulk auxiliary variables
    update_substep_auxils!(SOIL_BULK);

    # update the root auxiliary variables
    for root in ROOTS
        update_substep_auxils!(root);
    end;

    # update the junction auxiliary variables
    update_substep_auxils!(JUNCTION);

    # update the stem auxiliary variables
    update_substep_auxils!(TRUNK);
    for stem in BRANCHES
        update_substep_auxils!(stem);
    end;

    # update the leaf auxiliary variables
    for leaf in LEAVES
        update_substep_auxils!(leaf);
    end;

    return nothing
);

update_substep_auxils!(soil::SoilLayer{FT}) where {FT} = (
    soil.auxil.cp = heat_capacitance(soil);
    soil.auxil.t = soil.state.Σe / soil.auxil.cp;

    # update the conductance, potential, diffusivity, and thermal conductivity (0.5 for tortuosity factor)
    soil.auxil.k = relative_soil_k(soil.state.vc, soil.state.θ) * soil.state.vc.K_MAX * relative_viscosity(soil.auxil.t) / soil.auxil.δz;
    soil.auxil.ψ = soil_ψ_25(soil.state.vc, soil.state.θ; oversaturation = true) * relative_surface_tension(soil.auxil.t);
    soil.auxil.kd = 0.5 * max(0, soil.state.vc.Θ_SAT - soil.state.θ) / soil.auxil.δz;
    soil.auxil.kv = 0.5 * soil.state.vc.Θ_SAT / max(FT(0.01), soil.state.vc.Θ_SAT - soil.state.θ) / soil.auxil.δz;
    soil.auxil.λ_soil_water = (soil.state.λ_soil + soil.state.θ * Λ_THERMAL_H₂O(FT)) / soil.auxil.δz;

    # clean up the partial derivatives
    soil.auxil.∂e∂t = 0;
    soil.auxil.∂n∂t .= 0;
    soil.auxil.∂θ∂t = 0;

    return nothing
);

update_substep_auxils!(sbulk::SoilBulk{FT}) where {FT} = (
    # clear the dndt cahche
    sbulk.auxil.dndt .= 0;
    sbulk.auxil.runoff = 0;

    return nothing
);

update_substep_auxils!(root::Root{FT}) where {FT} = (
    # update root cp and temperature
    root.energy.auxil.cp = heat_capacitance(root);
    root.energy.auxil.t = root.energy.state.Σe / root.energy.auxil.cp;

    # update the root buffer pressure and flow
    x_aux = root.xylem.auxil;
    x_sta = root.xylem.state;
    if x_aux isa XylemHydraulicsAuxilNSS
        N = length(x_sta.v_storage);
        f_vis = relative_viscosity(root.energy.auxil.t);
        v_max_i = x_sta.v_max * x_sta.area * x_sta.l / N;
        for i in 1:N
            x_aux.p_storage[i] = capacitance_pressure(x_sta.pv, x_sta.v_storage[i] / v_max_i, root.energy.auxil.t);
            x_aux.flow_buffer[i] = (x_aux.p_storage[i] - (x_aux.pressure[i] + x_aux.pressure[i+1]) / 2) * x_sta.pv.k_refill / f_vis * x_sta.v_storage[i];
        end;
    end;

    # clear the partial derivatives
    root.energy.auxil.∂e∂t = 0;

    return nothing
);

update_substep_auxils!(junc::JunctionCapacitor{FT}) where {FT} = (
    # update junction cp and temperature
    junc.auxil.cp = heat_capacitance(junc);
    junc.auxil.t = junc.state.Σe / junc.auxil.cp;

    # update the junction buffer pressure
    junc.auxil.pressure = capacitance_pressure(junc.state.pv, junc.state.v_storage / junc.state.v_max, junc.auxil.t);

    # clear the partial derivatives
    junc.auxil.∂e∂t = 0;
    junc.auxil.∂w∂t = 0;

    return nothing
);

update_substep_auxils!(stem::Stem{FT}) where {FT} = (
    # update stem cp and temperature
    stem.energy.auxil.cp = heat_capacitance(stem);
    stem.energy.auxil.t = stem.energy.state.Σe / stem.energy.auxil.cp;

    # update the stem buffer pressure and flow
    x_aux = stem.xylem.auxil;
    x_sta = stem.xylem.state;
    if x_aux isa XylemHydraulicsAuxilNSS
        N = length(x_sta.v_storage);
        f_vis = relative_viscosity(stem.energy.auxil.t);
        v_max_i = x_sta.v_max * x_sta.area * x_sta.l / N;
        for i in 1:N
            x_aux.p_storage[i] = capacitance_pressure(x_sta.pv, x_sta.v_storage[i] / v_max_i, stem.energy.auxil.t);
            x_aux.flow_buffer[i] = (x_aux.p_storage[i] - (x_aux.pressure[i] + x_aux.pressure[i+1]) / 2) * x_sta.pv.k_refill / f_vis * x_sta.v_storage[i];
        end;
    end;

    # clear the partial derivatives
    stem.energy.auxil.∂e∂t = 0;

    return nothing
);

update_substep_auxils!(leaf::Leaf{FT}) where {FT} = (
    # update leaf cp and temperature
    leaf.energy.auxil.cp = heat_capacitance(leaf);
    leaf.energy.auxil.t = leaf.energy.state.Σe / leaf.energy.auxil.cp;

    # update the leaf buffer pressure and flow
    x_aux = leaf.xylem.auxil;
    c_aux = leaf.capacitor.auxil;
    if x_aux isa XylemHydraulicsAuxilNSS
        x_sta = leaf.xylem.state;
        c_sta = leaf.capacitor.state;
        f_vis = relative_viscosity(leaf.energy.auxil.t);
        c_aux.p = capacitance_pressure(c_sta.pv, c_sta.v_storage / c_sta.v_max, leaf.energy.auxil.t);
        c_aux.flow = (c_aux.p - x_aux.pressure[end]) * c_sta.pv.k_refill / f_vis * c_sta.v_storage * x_sta.area;
    else
        c_aux.flow = 0;
    end;

    # clear the partial derivatives
    leaf.energy.auxil.∂e∂t = 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_step_auxils! function
#
#######################################################################################################################################################################################################
"""

    update_step_auxils!(spac::MultiLayerSPAC{FT}) where {FT}

Update the auxiliary variables at big time step, given
- `spac` `MultiLayerSPAC` SPAC

"""
function update_step_auxils! end;

update_step_auxils!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; LEAVES) = spac;

    for leaf in LEAVES
        update_substep_auxils!(leaf);
    end;

    return nothing
);

update_step_auxils!(leaf::Leaf{FT}) where {FT} = (
    leaf.∫∂w∂t_out = 0;

    return nothing
);

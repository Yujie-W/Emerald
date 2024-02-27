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
#     2023-Oct-16: add leaf stomatal conductance variables calculations
#
#######################################################################################################################################################################################################
"""

    update_substep_auxils!(spac::BulkSPAC{FT}) where {FT}

Update the auxiliary variables at sub time step within a big time step, given
- `spac` `BulkSPAC` SPAC

"""
function update_substep_auxils! end;

update_substep_auxils!(spac::BulkSPAC{FT}) where {FT} = (
    soils = spac.soils;
    sbulk = spac.soil_bulk;

    roots = spac.plant.roots;
    junction = spac.plant.junction;
    trunk = spac.plant.trunk;
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;

    # update the soil layer auxiliary variables
    for soil in soils
        update_substep_auxils!(soil);
    end;

    # update soil bulk auxiliary variables
    update_substep_auxils!(sbulk);

    # update the root auxiliary variables
    for root in roots
        update_substep_auxils!(root);
    end;

    # update the junction auxiliary variables
    update_substep_auxils!(junction);

    # update the stem auxiliary variables
    update_substep_auxils!(trunk);
    for stem in branches
        update_substep_auxils!(stem);
    end;

    # update the leaf auxiliary variables
    for leaf in leaves
        update_substep_auxils!(leaf);
    end;

    return nothing
);

update_substep_auxils!(soil::SoilLayer{FT}) where {FT} = (
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
    # update the root buffer pressure and flow
    x_aux = root.xylem.auxil;
    x_tra = root.xylem.trait;
    x_sta = root.xylem.state;
    if x_aux isa XylemHydraulicsAuxilNSS
        N = length(x_sta.v_storage);
        f_vis = relative_viscosity(root.energy.s_aux.t);
        v_max_i = x_tra.v_max * x_tra.area * x_tra.l / N;
        for i in 1:N
            x_aux.p_storage[i] = capacitance_pressure(x_tra.pv, x_sta.v_storage[i] / v_max_i, root.energy.s_aux.t);
            x_aux.flow_buffer[i] = (x_aux.p_storage[i] - (x_aux.pressure[i] + x_aux.pressure[i+1]) / 2) * x_tra.pv.k_refill / f_vis * x_sta.v_storage[i];
        end;
    end;

    # clear the partial derivatives
    root.energy.auxil.∂e∂t = 0;

    return nothing
);

update_substep_auxils!(junc::JunctionCapacitor{FT}) where {FT} = (
    # clear the partial derivatives
    junc.auxil.∂e∂t = 0;
    junc.auxil.∂w∂t = 0;

    return nothing
);

update_substep_auxils!(stem::Stem{FT}) where {FT} = (
    # update the stem buffer pressure and flow
    x_aux = stem.xylem.auxil;
    x_tra = stem.xylem.trait;
    x_sta = stem.xylem.state;
    if x_aux isa XylemHydraulicsAuxilNSS
        N = length(x_sta.v_storage);
        f_vis = relative_viscosity(stem.energy.s_aux.t);
        v_max_i = x_tra.v_max * x_tra.area * x_tra.l / N;
        for i in 1:N
            x_aux.p_storage[i] = capacitance_pressure(x_tra.pv, x_sta.v_storage[i] / v_max_i, stem.energy.s_aux.t);
            x_aux.flow_buffer[i] = (x_aux.p_storage[i] - (x_aux.pressure[i] + x_aux.pressure[i+1]) / 2) * x_tra.pv.k_refill / f_vis * x_sta.v_storage[i];
        end;
    end;

    # clear the partial derivatives
    stem.energy.auxil.∂e∂t = 0;

    return nothing
);

update_substep_auxils!(leaf::Leaf{FT}) where {FT} = (
    # update the leaf buffer pressure and flow
    x_aux = leaf.xylem.auxil;
    c_aux = leaf.capacitor.auxil;
    if x_aux isa XylemHydraulicsAuxilNSS
        x_tra = leaf.xylem.trait;
        c_tra = leaf.capacitor.trait;
        c_sta = leaf.capacitor.state;
        f_vis = relative_viscosity(leaf.energy.s_aux.t);
        c_aux.p = capacitance_pressure(c_tra.pv, c_sta.v_storage / c_tra.v_max, leaf.energy.s_aux.t);
        c_aux.flow = (c_aux.p - x_aux.pressure[end]) * c_tra.pv.k_refill / f_vis * c_sta.v_storage * x_tra.area;
    else
        c_aux.flow = 0;
    end;

    # update the stomatal conductance
    leaf.flux.auxil.g_CO₂_shaded = 1 / (1 / leaf.flux.auxil.g_CO₂_b + FT(1.6) / leaf.flux.state.g_H₂O_s_shaded);
    for i in eachindex(leaf.flux.state.g_H₂O_s_sunlit)
        leaf.flux.auxil.g_CO₂_sunlit[i] = 1 / (1 / leaf.flux.auxil.g_CO₂_b + FT(1.6) / leaf.flux.state.g_H₂O_s_sunlit[i]);
    end;

    # clear the partial derivatives
    leaf.energy.auxil.∂e∂t = 0;
    leaf.flux.auxil.∂g∂t_shaded = 0;
    leaf.flux.auxil.∂g∂t_sunlit .= 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_step_auxils! function
#     2024-Jan-24: update leaf boundary layer conductance based on wind speed and leaf width
#
#######################################################################################################################################################################################################
"""

    update_step_auxils!(spac::BulkSPAC{FT}) where {FT}

Update the auxiliary variables at big time step, given
- `spac` `BulkSPAC` SPAC

"""
function update_step_auxils! end;

update_step_auxils!(spac::BulkSPAC{FT}) where {FT} = (
    airs = spac.airs;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;

    # clear the integrated values
    for leaf in leaves
        update_step_auxils!(leaf);
    end;

    # update leaf boundary layer conductance
    for i in eachindex(leaves)
        leaf = leaves[i];
        air = airs[lindex[i]];
        leaf.flux.auxil.g_CO₂_b = FT(0.14) * sqrt(air.auxil.wind / (FT(0.72) * leaf.bio.trait.width));
    end;

    return nothing
);

update_step_auxils!(leaf::Leaf{FT}) where {FT} = (
    leaf.flux.auxil.∫∂w∂t_out = 0;

    return nothing
);

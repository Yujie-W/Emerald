# This file contains functions to update the auxiliary variables at sub time steps

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_substep_auxils! function
#     2023-Sep-30: add buffer pressure and flow calculations
#     2023-Oct-06: add soil auxiliary variable calculations
#     2023-Oct-07: add soil bulk auxiliary variable calculations
#     2023-Oct-16: add leaf stomatal conductance variables calculations
#     2023-Oct-17: rename the function to substep_aux! to be consistent with the other function names such as t_aux!, s_aux!, dull_aux!, and step_aux!
#     2024-Jul-24: add leaf shedded flag
#     2024-Jul-30: compute OCS conductance along with CO₂ conductance
#
#######################################################################################################################################################################################################
"""

substep_aux!(spac::BulkSPAC{FT}) where {FT}

Update the auxiliary variables at sub time step within a big time step, given
- `spac` `BulkSPAC` SPAC

"""
function substep_aux! end;

substep_aux!(spac::BulkSPAC{FT}) where {FT} = (
    soils = spac.soils;
    sbulk = spac.soil_bulk;
    roots = spac.plant.roots;
    junction = spac.plant.junction;
    trunk = spac.plant.trunk;
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;

    # update the soil layer auxiliary variables
    for soil in soils
        substep_aux!(soil);
    end;

    # update soil bulk auxiliary variables
    substep_aux!(sbulk);

    # update the root auxiliary variables
    for root in roots
        substep_aux!(root);
    end;

    # update the junction auxiliary variables
    substep_aux!(junction);

    # update the stem auxiliary variables
    substep_aux!(trunk);
    for stem in branches
        substep_aux!(stem);
    end;

    # update the leaf auxiliary variables
    if !spac.plant._leaf_shedded
        for leaf in leaves
            substep_aux!(leaf);
        end;
    end;

    return nothing
);

substep_aux!(soil::SoilLayer{FT}) where {FT} = (
    # clean up the partial derivatives
    soil.auxil.∂e∂t = 0;
    soil.auxil.∂n∂t .= 0;
    soil.auxil.∂θ∂t = 0;

    return nothing
);

substep_aux!(sbulk::SoilBulk{FT}) where {FT} = (
    # clear the dndt cahche
    sbulk.auxil.dndt .= 0;
    sbulk.auxil.runoff = 0;

    return nothing
);

substep_aux!(root::Root{FT}) where {FT} = (
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

substep_aux!(junc::JunctionCapacitor{FT}) where {FT} = (
    # clear the partial derivatives
    junc.auxil.∂e∂t = 0;
    junc.auxil.∂w∂t = 0;

    return nothing
);

substep_aux!(stem::Stem{FT}) where {FT} = (
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

substep_aux!(leaf::CanopyLayer{FT}) where {FT} = (
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
    leaf.flux.auxil.g_CO₂ .= 1 ./ (1 ./ leaf.flux.auxil.g_CO₂_b .+ FT(1.6) ./ leaf.flux.state.g_H₂O_s);
    leaf.flux.auxil.g_OCS .= 1 ./ (1 ./ leaf.flux.auxil.g_OCS_b .+ FT(1.934) ./ leaf.flux.state.g_H₂O_s .+ 1 ./ (leaf.photosystem.trait.K_OCS .* leaf.photosystem.auxil.v_cmax));

    # clear the partial derivatives
    leaf.energy.auxil.∂e∂t = 0;
    leaf.energy.auxil.∂e∂t_le = 0;
    leaf.energy.auxil.∂e∂t_sh = 0;
    leaf.flux.auxil.∂g∂t .= 0;

    return nothing
);

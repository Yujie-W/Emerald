# This file contains functions to update the auxiliary variables at different time steps
#     - Big time step (such as the integrated values: GPP, ET, etc.)
#     - Sub time step within a big time step (such as the cache variables which can be recomputed: t, cp, and p)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add update_substep_auxils! function
#
#######################################################################################################################################################################################################
"""

    update_substep_auxils!(spac::MultiLayerSPAC{FT}) where {FT}

Update the auxiliary variables at sub time step within a big time step, given
- `spac` `MultiLayerSPAC` SPAC

"""
function update_substep_auxils! end;

update_substep_auxils!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; ROOTS, JUNCTION, TRUNK, BRANCHES, LEAVES) = spac;

    # update the soil auxiliary variables
    for soil in spac.SOIL.LAYERS
        update_substep_auxils!(soil);
    end;

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
    soil.∂e∂t = 0;
    soil.∂n∂t .= 0;
    soil.∂θ∂t = 0;

    return nothing
);

update_substep_auxils!(root::Root{FT}) where {FT} = (
    root.energy.auxil.cp = heat_capacitance(root);
    root.energy.auxil.t = root.energy.state.Σe / root.energy.auxil.cp;
    root.energy.auxil.∂e∂t = 0;

    return nothing
);

update_substep_auxils!(junc::JunctionCapacitor{FT}) where {FT} = (
    junc.auxil.cp = heat_capacitance(junc);
    junc.auxil.t = junc.state.Σe / junc.energy.auxil.cp;
    junc.auxil.∂e∂t = 0;
    junc.auxil.∂w∂t = 0;

    return nothing
);

update_substep_auxils!(stem::Stem{FT}) where {FT} = (
    stem.energy.auxil.cp = heat_capacitance(stem);
    stem.energy.auxil.t = stem.energy.state.Σe / stem.energy.auxil.cp;
    stem.energy.auxil.∂e∂t = 0;

    return nothing
);

update_substep_auxils!(leaf::Leaves2D{FT}) where {FT} = (
    leaf.NS.energy.auxil.cp = heat_capacitance(leaf);
    leaf.NS.energy.auxil.t = leaf.NS.energy.state.Σe / leaf.NS.energy.auxil.cp;
    leaf.NS.energy.auxil.∂e∂t = 0;

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

update_step_auxils!(leaf::Leaves2D{FT}) where {FT} = (
    leaf.∫∂w∂t_out = 0;

    return nothing
);

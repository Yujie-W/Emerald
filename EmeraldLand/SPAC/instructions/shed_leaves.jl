#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-24: add leaf shedding function
#     2024-Aug-06: set g_H₂O_s to 0 when leaves are shedded (so g will be 0 when regrowing)
#     2024-Sep-03: set asap to 0 when leaves are shedded
#     2024-Nov-05: move leaf shedding warning message into the function
#
#######################################################################################################################################################################################################
"""

    shed_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Shed leaves from the plant if ALLOW_LEAF_SHEDDING is true, given
- `config` Configuration struct
- `spac` `BulkSPAC` SPAC

"""
function shed_leaves! end;

shed_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    if !config.ALLOW_LEAF_SHEDDING || spac.plant._leaf_shedded
        return nothing
    end;

    # if leaf shedding is allowed, remember to update leaf area when roots are reconnected to the soil
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    can_str.trait.lai = 0;
    can_str.trait.δlai .= 0;
    for i in 1:n_layer
        leaf = leaves[i];
        leaf.xylem.trait.area = 0;
        leaf.xylem.state.asap = 0;
        leaf.flux.state.g_H₂O_s .= 0;
    end;

    # update the leaf shedding flag
    spac.plant._leaf_shedded = true;

    # update the canopy structure auxilary variables
    t_aux!(config, spac.canopy, spac.cache);

    # pop out warning message
    if config.MESSAGE_LEVEL == 2
        @warn "Leaf shedding is triggered";
    end;

    return nothing
);

#=
shed_leaves!(spac::BulkSPAC{FT}, lai_diff::FT) where {FT} = (
    @assert lai_diff < 0 "lai_diff should be negative when partially shedding leaves";

    can_str = spac.canopy.structure;
    junc = spac.plant.junction;
    leaves = spac.plant.leaves;
    sbulk = spac.bulk;
    n_layer = length(leaves);

    w_to_junc = 0;
    e_to_junc = 0;
    can_str.trait.lai += lai_diff;
    can_str.trait.δlai = can_str.trait.lai .* ones(FT, n_layer) ./ n_layer;
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        leaf = leaves[ilf];
        leaf.xylem.trait.area = sbulk.trait.area * can_str.trait.δlai[irt];
        w_to_junc -= sbulk.trait.area * lai_diff / n_layer * leaf.capacitor.state.v_storage;
        e_to_junc -= sbulk.trait.area * lai_diff / n_layer * leaf.capacitor.state.v_storage * CP_L_MOL(FT) * leaf.energy.s_aux.t;
    end;

    # update the junction state variables
    junc.state.v_storage += w_to_junc;
    junc.state.Σe += e_to_junc;

    # TODO: add the LMA to top soil

    return nothing
);
=#

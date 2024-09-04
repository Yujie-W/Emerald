#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-06: add function to regrow leaves (flag only)
#
#######################################################################################################################################################################################################
"""

    regrow_leaves_flag!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Regrow leaves if the following conditions are met, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` SPAC

"""
function regrow_leaves_flag!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    if !config.ALLOW_LEAF_SHEDDING || !spac.plant._leaf_shedded || spac.plant.pool.c_pool <= 0
        return nothing
    end;

    # set the regrow flag to true only if all roots are connected to the soil and the junction pressure is not too low
    flag = true;
    one_root = spac.plant.roots[1];
    p_50 = xylem_pressure(one_root.xylem.trait.vc, FT(0.5)) * relative_surface_tension(one_root.energy.s_aux.t);
    for s in spac.soils
        if s.s_aux.ψ < p_50
            flag = false;
            break;
        end;
    end;
    if (spac.plant.junction.s_aux.pressure < -0.1)
        flag = false;
    end;

    spac.plant._leaf_regrow = flag;

    return nothing
end;

#=
function grow_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, lai_diff::FT) where {FT}
    @assert lai_diff > 0 "lai_diff should be positive when regrowing leaves";

    can_str = spac.canopy.structure;
    junc = spac.plant.junction;
    leaves = spac.plant.leaves;
    sbulk = spac.bulk;
    n_layer = length(leaves);

    can_str.trait.lai += lai_diff;
    can_str.trait.δlai = can_str.trait.lai .* ones(FT, n_layer) ./ n_layer;
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        leaf = leaves[ilf];
        leaf.xylem.trait.area = sbulk.trait.area * can_str.trait.δlai[irt];
        delta_w = sbulk.trait.area * lai_diff / n_layer * leaf.capacitor.state.v_storage;
        delta_e = sbulk.trait.area * lai_diff / n_layer * leaf.capacitor.state.v_storage * CP_L_MOL(FT) * junc.energy.s_aux.t;
        junc.state.v_storage -= delta_w;
        junc.state.Σe -= delta_e;
    end;

    return nothing
end;
=#

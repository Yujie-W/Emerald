#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-06: add function to regrow leaves
#
#######################################################################################################################################################################################################
"""

    regrow_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Regrow leaves if the following conditions are met, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` SPAC

"""
function regrow_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    if !config.ALLOW_LEAF_SHEDDING || !spac.plant._leaf_shedded
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

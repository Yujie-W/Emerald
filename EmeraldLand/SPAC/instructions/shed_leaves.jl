#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-24: add leaf shedding function
#     2024-Aug-06: set g_H₂O_s to 0 when leaves are shedded (so g will be 0 when regrowing)
#
#######################################################################################################################################################################################################
"""

    shed_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Shed leaves from the plant if ALLOW_LEAF_SHEDDING is true, given
- `config` Configuration struct
- `spac` `BulkSPAC` SPAC

"""
function shed_leaves!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
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
        leaf.flux.state.g_H₂O_s .= 0;
    end;

    # update the canopy structure auxilary variables
    t_aux!(config, spac.canopy, spac.cache);

    return nothing
end;

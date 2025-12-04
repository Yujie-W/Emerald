#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-30: add function to push the mean temperature of the leaves to the history
#
#######################################################################################################################################################################################################
"""

    push_t_history!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT})

Push the mean temperature of the leaves to the history, given
- `config` SPAC configuration
- `spac` SPAC system

"""
function push_t_history!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # compute the mean temperature of the leaves
    sum_t::FT = 0;
    for l in spac.plant.leaves
        sum_t += l.energy.s_aux.t;
    end;
    mean_t = sum_t / length(spac.plant.leaves);

    # push the mean temperature to the history
    spac.plant.memory.t_history[spac.plant.memory.i_history] = mean_t;
    spac.plant.memory.i_history += 1;
    if spac.plant.memory.i_history > config.T_CLM_N
        spac.plant.memory.i_history = 1;
    end;

    return nothing
end;

# This file contains function to iterate the functions to make sure all time elapses per time step

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-18: add function time_stepper!
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2022-Nov-18: add option p_on to enable/disable plant flow and pressure profiles
#     2023-Jun-15: add judge for root connection
#     2023-Oct-18: redesign the logic flow to call the sub step functions in the time stepper
#     2024-Feb-29: when totoal count exceeds 1000, break the loop with an error (this need to be set as a BUG instead of warning)
#     2024-Aug-06: move leaf shedding condition into the time_stepper function (otherwise leaf shedding will be triggered immediately after the regrowth because the pressure is not reset)
#     2024-Nov-05: move leaf shedding warning message into the function call (to avoid the warning message when the leaf is not shedded)
#
#######################################################################################################################################################################################################
"""

    time_stepper!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::Number) where {FT}

Move forward in time for SPAC with time stepper controller, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC
- `δt` Time step (if not given, solve for steady state solution)

"""
function time_stepper!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::Number) where {FT}
    # run the update function until time elapses
    count = 0;
    δt_remain::FT = δt;
    while δt_remain > 0
        substep_preparations!(config, spac);
        δt_step = adjusted_time(spac, δt_remain);
        substep_budgets!(config, spac, δt_step);
        δt_remain -= δt_step;

        # determine whether to shed leaves at the end of each sub time step
        bottom_leaf = spac.plant.leaves[1];
        p_crt = xylem_pressure(bottom_leaf.xylem.trait.vc, config.KR_THRESHOLD) * relative_surface_tension(bottom_leaf.energy.s_aux.t);
        if !spac.plant._leaf_shedded && bottom_leaf.xylem.auxil.pressure[end] < p_crt
            shed_leaves!(config, spac);
        end;

        # if total count exceeds 1000, break the loop
        count += 1;
        if (count > 1000) && (δt_step < 0.01) && (δt_remain > 10)
            @error "Number of steppers exceeds 1000, breaking..." δt_remain δt_step;
            return error("Number of steppers exceeds 1000!")
        end;
    end;

    return nothing
end;

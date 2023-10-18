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
    t_res = FT(δt);
    while true
        count += 1;
        t_step = adjusted_time(config, spac, t_res);

        # run the budgets for all ∂x∂t
        soil_budgets!(config, spac, t_step);
        stomatal_conductance!(spac, t_step);
        spac_energy_budget!(spac, t_step);
        if spac.plant._root_connection
            plant_water_budget!(spac, t_step);
        end;

        t_res -= t_step;

        # if _t_res > 0 rerun the budget functions (shortwave radiation not included) and etc., else break
        if t_res > 0
            update_substep_auxils!(spac);
            longwave_radiation!(config, spac);
            if spac.plant._root_connection
                plant_flow_profile!(config, spac);
                plant_pressure_profile!(config, spac);
            end;
            plant_photosynthesis!(spac, GCO₂Mode());
            soil_profiles!(config, spac);
            stomatal_conductance_profile!(spac);
            spac_energy_flow!(spac);
        else
            break;
        end;

        # if total count exceeds 100
        if (count > 1000) && (t_step < 0.01) && (t_res > 10)
            @info "Number of steppers exceeds 1000, breaking..." spac.canopy.structure.state.lai t_res t_step;
            break;
        end;
    end;

    return nothing
end;

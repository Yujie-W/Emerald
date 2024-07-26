# this file contains functions to run at sub steps

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Feb-27: add function substep_preparations! to run the steps to compute partial deriative of the state variables
#     2024-Feb-28: add β_factor! to compute β factor before running plant_photosynthesis! and stomatal_conductance_profile! (for empirical stomatal models)
#
#######################################################################################################################################################################################################
"""

    substep_preparations!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Run the steps to compute partial deriative of the state variables, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC

"""
function substep_preparations!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    substep_aux!(spac);
    longwave_radiation!(spac);
    plant_flow_profile!(config, spac);
    plant_pressure_profile!(config, spac);      # 3rd most time consuming function (about 14%)
    soil_profiles!(config, spac);
    β_factor!(spac);
    @time plant_photosynthesis!(spac, GCO₂Mode());    # the most time consuming function (about 59%) # TODO
    stomatal_conductance_profile!(spac);        # 2nd most time consuming function (about 27%)
    spac_energy_flow!(config, spac);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Feb-27: add function substep_budgets! to run the budgets for all ∂x∂t using the adjusted time step
#
#######################################################################################################################################################################################################
"""

    substep_budgets!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::FT) where {FT}

Run the budgets for all ∂x∂t using the adjusted time step, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC
- `δt` Time step

"""
function substep_budgets!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::FT) where {FT}
    soil_budgets!(config, spac, δt);
    plant_water_budget!(spac, δt);
    spac_energy_budget!(spac, δt);
    stomatal_conductance!(spac, δt);
    s_aux!(spac);

    return nothing
end;

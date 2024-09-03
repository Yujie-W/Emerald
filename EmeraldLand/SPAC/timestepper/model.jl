# This file contains functions to run the SPAC model in a large time step (update_step_auxils! will be called)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to run soil plant air continuum at a time step
#     2022-Jul-13: add soil water and energy budget
#     2022-Aug-11: run fluorescence model after photosynthesis
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2022-Nov-18: add option p_on to enable/disable plant flow and pressure profiles
#     2023-Mar-28: if no root is connected, set LAI = 0
#     2023-Mar-28: run PlantHydraulics as the first step
#     2023-Apr-08: set runoff to 0 at the beginning of each time interval
#     2023-Jun-15: add judge for root connection (avoid a bug in non-vegetated land)
#     2023-Oct-18: design the logic flow with new time_stepper! function design
#     2024-Jul-24: add leaf shedded flag
#     2024-Jul-24: set regrow threshold to 50% of the critical pressure
#     2024-Aug-06: move leaf shedding condition into the substep_aux! function
#     2024-Aug-06: set regrow threshold to when junction pressure is higher than -0.1 MPa
#     2024-Sep-03: add step to grow xylem when carbon pool is higher than the threshold
# To do
#     TODO: add top soil evaporation
#
#######################################################################################################################################################################################################
"""
This function runs the model using the following steps:
- Run hydraulic model, and determine whether to set LAI = 0
- Run canopy RT model
- Run photosynthesis model
- Run canopy fluorescence model
- Run soil water and energy budget (calculate ∂Θ∂t and ∂e∂t only)
- Run leaf stomatal conductances (calculate ∂g∂t only)
- Run leaf energy budget (calculate ∂T∂t only)
- Run time stepper (using ∂X∂t * δt, and make sure δt is not too high)

This function is supposed to have the highest hierarchy, and should support all SPAC types defined in EmeraldNamespace.jl. Note to update the water flow profile when initializing the SPAC.

    soil_plant_air_continuum!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::Number) where {FT}

Run SPAC model and move forward in time with time stepper controller, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC,
- `δt` Time step

"""
function soil_plant_air_continuum!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::Number) where {FT}
    # 1. run the functions are do not need to be run at sub time step (e.g. shortwave radiation)
    step_preparations!(config, spac);

    # 2. use time stepper to run the functions that need to be run at sub time step
    #    this function takes almost all the time, go deep to see what make it so slow
    time_stepper!(config, spac, δt);

    # 3. run canopy reflectance and fluorescence to use with remote sensing
    step_remote_sensing!(config, spac);

    # 4. update the legacy for the next time step
    update_legacy!(config, spac);

    # 5. determine whether to regrow the leaves in the next round of LAI update (set flag only)
    regrow_leaves_flag!(config, spac);

    # 6. determine whether tree death occurs
    plant_death!(config, spac);

    # 7. grow xylem when carbon pool is higher than the threshold
    plant_growth!(spac);

    return nothing
end;

# add an alias for soil_plant_air_continuum!
"""

    spac!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::Number) where {FT}

Run SPAC model and move forward in time with time stepper controller, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC,
- `δt` Time step

"""
spac! = soil_plant_air_continuum!;

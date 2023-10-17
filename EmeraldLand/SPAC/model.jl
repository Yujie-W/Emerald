#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to run soil plant air continuum at a time step
#     2022-Jul-13: use soil_energy! and soil_water! for soil water and energy budget
#     2022-Aug-11: run fluorescence model after photosynthesis
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Sep-07: add method to solve for steady state solution
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2022-Nov-18: add option p_on to enable/disable plant flow and pressure profiles
#     2023-Mar-11: add method for a SPAC == nothing
#     2023-Mar-28: if no root is connected, set LAI = 0
#     2023-Mar-28: run PlantHydraulics as the first step
#     2023-Apr-08: set runoff to 0 at the beginning of each time interval
#     2023-Apr-13: add config to function call
#     2023-Jun-15: add judge for root connection (avoid a bug in non-vegetated land)
#     2023-Sep-07: add ALLOW_LEAF_SHEDDING check
#     2023-Sep-11: move the optional t_on and θ_on to the config struct
#     2023-Sep-14: remove some if else control from root disconnection
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

    soil_plant_air_continuum!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::Number) where {FT}
    soil_plant_air_continuum!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    soil_plant_air_continuum!(config::Nothing, spac::Nothing, δt::Number) where {FT}
    soil_plant_air_continuum!(config::Nothing, spac::Nothing) where {FT}

Run SPAC model and move forward in time with time stepper controller, given
- `spac` `MultiLayerSPAC` SPAC, or nothing
- `config` Configuration for `MultiLayerSPAC`
- `δt` Time step (if not given, solve for steady state solution)

"""
function soil_plant_air_continuum! end;

# TODO: add lite mode later to update energy balance (only longwave radiation and soil+leaf energy budgets)? Or use shorter time steps (will be time consuming, but more accurate)
# TODO: add top soil evaporation
soil_plant_air_continuum!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::Number) where {FT} = (
    # 0. set total runoff to 0 so as to accumulate with sub-timestep
    spac.SOIL_BULK.auxil.runoff = 0;
    update_substep_auxils!(spac);

    # 1. run plant hydraulic model (must be run before plant_photosynthesis! as the latter may need β for empirical models)
    plant_flow_profile!(config, spac);
    plant_pressure_profile!(config, spac);
    (!spac._root_connection && config.ALLOW_LEAF_SHEDDING) ? update!(config, spac; lai = 0) : nothing;

    # 2. run canopy RT
    canopy_radiation!(config, spac);

    # 3. run photosynthesis model
    plant_photosynthesis!(spac, GCO₂Mode());

    # save the result at this stage for the results at the beginning of this time step

    # 5. run soil energy water budget
    soil_profiles!(config, spac);

    # 6. run leaf stomatal conductance budget
    stomatal_conductance_profile!(spac);

    # 7. run plant energy budget
    spac_energy_flow!(spac);

    # 8. update the prognostic variables
    time_stepper!(config, spac, δt);

    # 4. run canopy reflectance and fluorescence
    sensor_geometry!(config, spac);
    reflection_spectrum!(config, spac);
    fluorescence_spectrum!(config, spac);

    return nothing
);

soil_plant_air_continuum!(config::Nothing, spac::Nothing, δt::Number) = nothing;


# add an alias for soil_plant_air_continuum!
spac! = soil_plant_air_continuum!;

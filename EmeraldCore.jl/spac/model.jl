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
#     2023-Apr-13: add config to function call
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

    soil_plant_air_continuum!(
                spac::MultiLayerSPAC{FT},
                config::MultiLayerSPACConfiguration{FT},
                δt::Number;
                p_on::Bool = true,
                t_on::Bool = true,
                update::Bool = false,
                θ_on::Bool = true) where {FT<:AbstractFloat}
    soil_plant_air_continuum!(spac::MultiLayerSPAC{FT}, config::MultiLayerSPACConfiguration{FT}; update::Bool = false) where {FT<:AbstractFloat}
    soil_plant_air_continuum!(spac::Nothing, config::Nothing, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat}
    soil_plant_air_continuum!(spac::Nothing, config::Nothing; update::Bool = false) where {FT<:AbstractFloat}

Run SPAC model and move forward in time with time stepper controller, given
- `spac` `MultiLayerSPAC` SPAC, or nothing
- `config` Configuration for `MultiLayerSPAC`
- `δt` Time step (if not given, solve for steady state solution)
- `p_on` If true, plant hydraulic flow and pressure profiles will be updated
- `t_on` If true, plant energy budget is on (set false to run sensitivity analysis or prescribing mode)
- `θ_on` If true, soil water budget is on (set false to run sensitivity analysis or prescribing mode)
- `update` If true, update leaf xylem legacy effect

"""
function soil_plant_air_continuum! end

# TODO: add lite mode later to update energy balance (only longwave radiation and soil+leaf energy budgets)? Or use shorter time steps (will be time consuming, but more accurate)
# TODO: add top soil evaporation
soil_plant_air_continuum!(spac::MultiLayerSPAC{FT}, config::MultiLayerSPACConfiguration{FT}, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat} = (
    # 1. run plant hydraulic model (must be run before leaf_photosynthesis! as the latter may need β for empirical models)
    xylem_flow_profile!(spac, FT(0));
    xylem_pressure_profile!(spac; update = update);
    if !spac._root_connection
        update!(spac; lai = 0);

        return nothing
    end;

    # 2. run canopy RT
    canopy_radiation!(spac);

    # 3. run photosynthesis model
    stomatal_conductance_profile!(spac);
    leaf_photosynthesis!(spac, GCO₂Mode());

    # 4. run canopy fluorescence
    canopy_fluorescence!(spac, config);

    # save the result at this stage for the results at the beginning of this time step

    # 5. run soil energy water budget
    soil_budget!(spac);

    # 6. run leaf stomatal conductance budget
    stomatal_conductance!(spac);

    # 7. run plant energy budget
    plant_energy!(spac);

    # 8. update the prognostic variables
    time_stepper!(spac, δt; p_on = p_on, t_on = t_on, update = update, θ_on = θ_on);

    return nothing
);

soil_plant_air_continuum!(spac::Nothing, config::Nothing, δt::Number; p_on::Bool = true, t_on::Bool = true, update::Bool = false, θ_on::Bool = true) = nothing;

soil_plant_air_continuum!(spac::MultiLayerSPAC{FT}, config::MultiLayerSPACConfiguration{FT}; update::Bool = false) where {FT<:AbstractFloat} = (
    # 1. run plant hydraulic model (must be run before leaf_photosynthesis! as the latter may need β for empirical models)
    xylem_flow_profile!(spac, FT(0));
    xylem_pressure_profile!(spac; update = update);
    if !spac._root_connection
        update!(spac; lai = 0);

        return nothing
    end;

    # 2. run canopy RT
    canopy_radiation!(spac);

    # 3. update the prognostic variables (except for soil water and temperature)
    time_stepper!(spac, config; update = update);

    # save the result at this stage for the results at the steady state

    return nothing
);

soil_plant_air_continuum!(spac::Nothing, config::Nothing; update::Bool = false) = nothing;


# add an alias for soil_plant_air_continuum!
spac! = soil_plant_air_continuum!;

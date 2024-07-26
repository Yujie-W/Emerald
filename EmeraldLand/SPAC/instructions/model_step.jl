# this file contains functions to run at a big step that will not be run at sub time steps to save time (e.g., shortwave radiation)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function to run preparation functions at the beginning of each time step
#     2024-Jun-06: compute the canopy structure every time (because soil moisture may change)
#
#######################################################################################################################################################################################################
"""

    step_preparations!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Run the functions that do not need to be run at sub time step (these functions are not impacted by states that change in substeps), given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC

"""
function step_preparations!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    step_aux!(spac);
    soil_albedo!(config, spac);
    canopy_structure!(config, spac);
    sun_geometry_aux!(config, spac);
    sun_geometry!(config, spac);
    shortwave_radiation!(config, spac);     # memory allocation due to broadcasting 2D to 1D
    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function to perform remote sensing analysis at the end of each time step
#
#######################################################################################################################################################################################################
"""

    step_remote_sensing!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Run the remote sensing functions that do not need to be run at sub time step (these functions may be impacted by states that change in substeps), given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC

"""
function step_remote_sensing!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    sensor_geometry_aux!(config, spac);
    sensor_geometry!(config, spac);
    reflection_spectrum!(config, spac);
    fluorescence_spectrum!(config, spac);

    return nothing
end;

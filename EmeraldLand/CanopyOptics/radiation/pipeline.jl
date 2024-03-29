# This file contains a pipeline of functions to compute the canopy radiation

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-14: redo canopy_radiation!
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update the canopy radiation related auxiliary variables, given
- `config` SPAC configuration
- `spac` SPAC

"""
function canopy_radiation!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # update the canopy structure
    inclination_angles!(config, spac);
    canopy_structure!(config, spac);
    soil_albedo!(config, spac);
    sun_geometry!(config, spac);

    # run longwave and shortwave radiation
    longwave_radiation!(config, spac);
    shortwave_radiation!(config, spac);

    return nothing
end;

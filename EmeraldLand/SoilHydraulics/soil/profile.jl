# This file contains the function to update the flow, diffusion, and energy flow profiles in the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-07: add function to update the flow, diffusion, and energy flow profiles in the soil
#
#######################################################################################################################################################################################################
"""

    soil_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update the flow, diffusion, and energy flow profiles in the soil, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` type SPAC

"""
function soil_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    soil_water_infiltration!(spac);
    root_source_sink!(spac);
    trace_gas_diffusion!(config, spac);

    return nothing
end;

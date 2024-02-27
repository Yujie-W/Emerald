# this file contains functions to update the slowly changing auxiliary variables that do not change at different time steps (to speed things up)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function to update the slowly changing auxiliary variables
#
#######################################################################################################################################################################################################
"""

    dull_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update the slowly changing auxiliary variables, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC

"""
function dull_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # the dull auxiliary variables include
    #     1 . leaf optical properties related to leaf pigments
    #     2 . canopy optical properties related to leaf pigments and canopy structure (e.g., canopy diffuse radiation related parameters)
    plant_leaf_spectra!(config, spac);
    canopy_structure!(config, spac);

    return nothing
end;

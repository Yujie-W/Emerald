# This file contains function to initialize the structs used in the SPAC model

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-27: add function to initialize SPAC
#     2022-Jun-27: add leaf area controller to make sure soil and leaf areas are consistent with leaf area index
#     2023-Mar-27: initialize soil and leaf e as well (because T, SWC may be changed)
#     2023-Jun-12: initialize soil trace gas as well
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-09: add root and stem initialization in the initialization of SPAC
#     2023-Oct-17: update step and subtep auxils during initialization
#     2023-Oct-18: initialize leaf inclination angles and canopy structure during initialization
#     2024-Feb-23: separate initialize_states! from initialize_spac!
#
#######################################################################################################################################################################################################
"""

    initialize_spac!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Initialize the SPAC, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` SPAC

"""
function initialize_spac!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    # make sure tha auxilary variables are correctly initialized
    update_step_auxils!(spac);
    update_substep_auxils!(spac);

    # initialize leaf level spectra
    plant_leaf_spectra!(config, spac);

    # initialize canopy structure
    t_aux!(config, spac.canopy.structure.trait, spac.canopy.structure.t_aux);
    canopy_structure!(config, spac);

    return nothing
end;

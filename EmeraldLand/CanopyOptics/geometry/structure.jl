# This file contains function to compute the canopy structure related parameters before running the canopy optical models

#######################################################################################################################################################################################################
#
# Changes to this file
# General
#     2023-Oct-10: add function canopy_structure! (run only once per inclination angle distribution)
#
#######################################################################################################################################################################################################
"""

    canopy_structure!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}

Update canopy structure related auxiliary variables, given
- `config` SPAC configuration
- `can` SPAC canopy

"""
function canopy_structure!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}
    (; Θ_INCL) = config;

    # compute the weighed average of the leaf inclination angle distribution
    can.structure.auxil.bf = 0;
    for i in eachindex(Θ_INCL)
        can.structure.auxil.bf += can.structure.state.p_incl[i] * cosd(Θ_INCL[i]) ^ 2;;
    end;

    return nothing
end;

# This file contains functions to compute leaf photosynthesis of the entire plant

######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for MultiLayerSPAC
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Mar-11: only compute respiration rate if solar zenith angle >= 89
#     2023-Mar-11: do nothing if LAI == 0
#
#######################################################################################################################################################################################################
"""

    plant_photosynthesis!(spac::MultiLayerSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `MultiLayerSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`

"""
function plant_photosynthesis!(spac::MultiLayerSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}
    (; AIRS, CANOPY, LEAVES, LEAVES_INDEX) = spac;

    if CANOPY.structure.state.lai == 0
        return nothing
    end;

    rd_only = CANOPY.sun_geometry.state.sza < 89 ? false : true;
    for i in eachindex(LEAVES)
        leaf_photosynthesis!(LEAVES[i], AIRS[LEAVES_INDEX[i]], mode; rd_only = rd_only);
    end;

    return nothing
end;

# This file contains functions to compute leaf photosynthesis of the entire plant

######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for BulkSPAC
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Mar-11: only compute respiration rate if solar zenith angle >= 89
#     2023-Mar-11: do nothing if LAI == 0
#
#######################################################################################################################################################################################################
"""

    plant_photosynthesis!(spac::BulkSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `BulkSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`

"""
function plant_photosynthesis!(spac::BulkSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}
    if spac.canopy.structure.state.lai == 0
        return nothing
    end;

    airs = spac.airs;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;

    rd_only = spac.canopy.sun_geometry.state.sza > 89;
    for i in eachindex(leaves)
        leaf_photosynthesis!(leaves[i], airs[lindex[i]], mode; rd_only = rd_only);
    end;

    return nothing
end;
